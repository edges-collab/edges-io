"""
This module defines the overall file structure and internal contents of the
calibration observations. It does *not* implement any algorithms/methods on that data,
making it easier to separate the algorithms from the data checking/reading.
"""
from __future__ import annotations

import attr
import logging
import numpy as np
import re
import read_acq
import tempfile
import toml
import warnings
import yaml
from bidict import bidict
from cached_property import cached_property
from copy import copy
from datetime import datetime
from io import StringIO
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

from . import utils
from ._structure import _DataContainer, _DataFile
from .data import DATA_PATH
from .h5 import HDF5RawSpectrum
from .logging import logger

with open(DATA_PATH / "calibration_loads.toml") as fl:
    data = toml.load(fl)
    LOAD_ALIASES = bidict({v["alias"]: k for k, v in data.items()})
    LOAD_MAPPINGS = {
        v: k
        for k, val in data.items()
        for v in val.get("misspells", []) + [val["alias"]]
    }

with open(DATA_PATH / "antenna_simulators.toml") as fl:
    ANTENNA_SIMULATORS = toml.load(fl)

# Dictionary of misspelled:true mappings.
ANTSIM_REVERSE = {
    v: k for k, val in ANTENNA_SIMULATORS.items() for v in val.get("misspells", [])
}


@attr.s
class _SpectrumOrResistance(_DataFile):
    load_pattern = "|".join(LOAD_ALIASES.values())
    antsim_pattern = "|".join(ANTENNA_SIMULATORS.keys())
    _antsim_rev_pattern = "|".join(ANTSIM_REVERSE.keys())
    _load_rev_pattern = "|".join(LOAD_MAPPINGS.keys())
    _loadname_pattern = (
        f"{load_pattern}|{antsim_pattern}|{_antsim_rev_pattern}|{_load_rev_pattern}"
    )

    pattern = (
        fr"(?P<load_name>{load_pattern}|{antsim_pattern})"
        + r"_(?P<run_num>\d{2})_(?P<year>\d{4})_(?P<day>\d{3})_("
        r"?P<hour>\d{2})_(?P<minute>\d{2})_(?P<second>\d{2})_lab.(?P<file_format>\w{2,"
        r"3})$"
    )
    write_pattern = (
        "{load_name}_{run_num:0>2}_{year:0>4}_{jd:0>3}_{hour:0>2}_{minute:0>2}_"
        "{second:0>2}_lab.{file_format}"
    )

    known_patterns = (
        (
            r"^(?P<load_name>%s)" % _loadname_pattern
            + r"_25C_(?P<month>\d{1,2})_(?P<day>\d{1,2})_("
            r"?P<year>\d\d\d\d)_(?P<hour>\d{1,2})_(?P<minute>\d{"
            r"1,2})_(?P<second>\d{1,2}).(?P<file_format>\w{2,3})$"
        ),
        (
            r"^(?P<load_name>%s)" % _loadname_pattern
            + r"_(?P<month>\d{1,2})_(?P<day>\d{1,2})_("
            r"?P<year>\d\d\d\d)_(?P<hour>\d{1,2})_(?P<minute>\d{"
            r"1,2})_(?P<second>\d{1,2}).(?P<file_format>\w{2,3})$"
        ),
        (
            f"^(?P<load_name>{_loadname_pattern})"
            + r"(?P<run_num>\d{1,2})_25C_(?P<month>\d{1,"
            r"2})_(?P<day>\d{1,2})_(?P<year>\d\d\d\d)_("
            r"?P<hour>\d{1,2})_(?P<minute>\d{1,"
            r"2})_(?P<second>\d{1,2}).(?P<file_format>\w{2,3})$"
        ),
        (
            r"(?P<load_name>%s)" % _loadname_pattern
            + r"_(?P<year>\d{4})_(?P<day>\d{3})_"
            r"(?P<hour>\d{2})_(?P<minute>\d{2})_(?P<second>\d{2})_lab."
            r"(?P<file_format>\w{2,3})$"
        ),
        (
            r"(?P<load_name>%s)" % _loadname_pattern
            + r"_(?P<year>\d{4})_(?P<day>\d{3})_(?P<hour>\d{2}).(?P<file_format>\w{2,3})$"
        ),
        (
            r"(?P<load_name>%s)" % _loadname_pattern
            + r"_(?P<year>\d{4})_(?P<day>\d{3})_(?P<hour>\d{2}).(?P<file_format>\w{2,3})$"
        ),
        (
            r"(?P<load_name>%s)" % _loadname_pattern
            + r"_(?P<year>\d{4})_(?P<day>\d{3})_lab.(?P<file_format>\w{2,3})$"
        ),
        (
            r"(?P<load_name>%s)" % _loadname_pattern
            + r"_(?P<run_num>\d)_(?P<year>\d{4})_(?P<day>\d{3})_lab.(?P<file_format>\w{2,"
            r"3})$"
        ),
        (
            r"(?P<load_name>%s)" % _loadname_pattern
            + r"_\d{2}C_(?P<month>\d{1,2})_(?P<day>\d{1,2})_(?P<year>\d{4})_(?P<hour>\d{"
            r"1,2})_(?P<minute>\d{1,2})_(?P<second>\d{1,2}).(?P<file_format>\w{2,3})"
        ),
        (
            r"(?P<load_name>%s)" % _loadname_pattern
            + r"_(?P<run_num>\d{2})_(?P<year>\d{4})_(?P<day>\d{3})_("
            r"?P<hour>\d{2})_(?P<minute>\d{2}).(?P<file_format>\w{2,3})$"
        ),
        (
            r"(?P<load_name>%s)" % _loadname_pattern
            + r"_(?P<run_num>\d{2})_(?P<year>\d{4})_(?P<day>\d{3})_("
            r"?P<hour>\d{2}).(?P<file_format>\w{2,3})$"
        ),
        (
            r"(?P<load_name>%s)" % _loadname_pattern
            + r"_(?P<year>\d{4})_(?P<day>\d{3})_("
            r"?P<hour>\d{2})_(?P<minute>\d{2}).(?P<file_format>\w{2,3})$"
        ),
        (
            r"(?P<load_name>%s)" % _loadname_pattern
            + r"_(?P<year>\d{4})_(?P<day>\d{3})_(?P<hour>\d{2}).(?P<file_format>\w{2,3})$"
        ),
    )

    known_substitutions = [
        ("degC", "C"),
        ("_25C", ""),
        ("_15C", ""),
        ("_35C", ""),
        ("LongCableShort_", "LongCableShorted_"),
    ]

    supported_formats = []

    @classmethod
    def typestr(cls, name: str):
        return cls.__name__ + re.match(cls.pattern, name).groupdict()["load_name"]

    @classmethod
    def _get_filename_parameters(cls, dct: dict):
        out = {"run_num": 1, "hour": 0, "minute": 0, "second": 0}

        if "month" in dct:
            out["jd"] = utils.ymd_to_jd(dct["year"], dct["month"], dct["day"])
        elif "day" in dct:
            out["jd"] = dct["day"]

        # Switch Antenna Simulator "misspells" to true form.
        if dct["load_name"] in ANTSIM_REVERSE:
            dct["load_name"] = ANTSIM_REVERSE[dct["load_name"]]

        elif dct["load_name"] in LOAD_MAPPINGS:
            dct["load_name"] = LOAD_MAPPINGS[dct["load_name"]]

        return out

    @classmethod
    def _validate_match(cls, match: dict[str, str], filename: str):
        if int(match["run_num"]) < 1:
            logger.error(f"The run_num for {filename} is less than one!")
        if not (2010 <= int(match["year"]) <= 2030):
            logger.error(f"The year for {filename} ({match['year']}) is a bit strange!")
        if not (0 <= int(match["day"]) <= 366):
            logger.error(
                f"The day for {filename} ({match['day']}) is outside the number "
                f"of days in a year"
            )
        if not (0 <= int(match["hour"]) <= 24):
            logger.error(f"The hour for {filename} is outside 0-24!")
        if not (0 <= int(match["minute"]) <= 60):
            logger.error(f"The minute for {filename} is outside 0-60!")
        if not (0 <= int(match["second"]) <= 60):
            logger.error(f"The second for {filename} is outside 0-60!")
        if match["file_format"] not in cls.supported_formats:
            logger.error(
                f"The file {filename} is not of a supported format "
                f"({cls.supported_formats}). Got format {match['file_format']}"
            )

    @classmethod
    def from_load(
        cls,
        load: str,
        direc: str | Path,
        run_num: int | None = None,
        filetype: str | None = None,
    ) -> list[_SpectrumOrResistance]:
        """
        Initialize the object in a simple way.

        Parameters
        ----------
        load
            The load name (eg. 'Ambient', 'HotLoad') or its alias (eg. 'ambient',
            'hot_load').
        direc
            The directory in which to search for relevant data
        run_num
            The run number of the data to use. Default, the last run. Each run is
            independent and different run_nums may be used for different loads.
        filetype
            The filetype of the data. Must be one of the supported formats. Defaults
            to `_default_filetype`.
        """
        direc = Path(direc)

        if load in LOAD_ALIASES:
            load = LOAD_ALIASES[load]

        if load not in LOAD_ALIASES.values() and load not in ANTENNA_SIMULATORS:
            logger.error(
                f"The load specified [{load}] is not one of the options available."
            )

        files = sorted(direc.glob(f"{load}_??_????_???_??_??_??_lab.*"))

        if not files:
            raise utils.LoadExistError(
                f"No files exist for the load {load} for any filetype on that path: {direc}."
                f"Found files: {list(files)}."
            )

        filetype = [filetype] if filetype else cls.supported_formats

        # Use any format so long as it is supported
        rfiles = []
        for ftype in filetype:
            rfiles = [fl for fl in files if fl.suffix == ("." + ftype)]
            if rfiles:
                break

        if not rfiles:
            raise utils.LoadExistError(
                f"No files exist for the load {load} for any of the filetypes '{filetype}'."
                f"Found files: {list(files)}."
            )

        files = rfiles

        # Restrict to the given run_num (default last run)
        run_nums = [int(fl.name[len(load) + 1 : len(load) + 3]) for fl in files]
        if run_num is None:
            run_num = max(run_nums)

        pre_files = files.copy()
        files = [fl for fl, num in zip(files, run_nums) if num == run_num]

        if not files:
            raise ValueError(
                f"No {load} files exist on path ({direc}) with run_num={run_num}. "
                f"Potential files: {pre_files}"
            )

        return [cls(fl) for fl in files]

    @cached_property
    def run_num(self):
        """The run number of the data. All run_nums must be the same for all files in
        the data.

        Every observation may have several runs. Note that different runs may be mixed
        for different loads.
        """
        # Ensure all load names are the same
        return self._match_dict["run_num"]

    @cached_property
    def year(self):
        """Year on which data acquisition began"""
        # Ensure all load names are the same
        return int(self._match_dict["year"])

    @cached_property
    def day(self) -> int:
        return int(self._match_dict["day"])

    @cached_property
    def load_name(self):
        return LOAD_ALIASES.inverse.get(
            self._match_dict["load_name"], self._match_dict["load_name"]
        )

    @cached_property
    def hour(self):
        """List of integer hours (one per file) at which data acquisition was begun"""
        return int(self._match_dict["hour"])

    @cached_property
    def minute(self):
        """List of integer minutes (one per file) at which data acquisition was begun"""
        return int(self._match_dict["minute"])

    @cached_property
    def second(self):
        """List of integer seconds (one per file) at which data acquisition was begun"""
        return int(self._match_dict["second"])


@attr.s
class FieldSpectrum:
    path: str | Path = attr.ib(converter=Path)

    @path.validator
    def _pth_vld(self, att, val):
        if not val.exists():
            raise OSError(f"{self.path} does not exist!")

        if self.file_format not in ["h5", "acq"]:
            raise TypeError(f"{self.path} has bad file format, must be h5 or acq")

    @cached_property
    def file_format(self) -> str:
        """The file format of the data to be read."""
        return self.path.suffix[1:]

    @cached_property
    def data(self) -> HDF5RawSpectrum | list[HDF5RawSpectrum]:
        """A view of the data in the file as a HDF5Object.

        If the file is an ACQ file, it will be read completely into memory and cast
        into the same format as a :class:`~h5.HDF5Object` so that the API is the same.

        If the number of files is more than one, `data` will be a list of objects.
        """
        if self.file_format == "h5":
            return HDF5RawSpectrum(self.path)
        elif self.file_format == "acq":
            spec, freq, time, meta = self._read_acq(self.path)
            return HDF5RawSpectrum.from_data(
                {
                    "spectra": spec,
                    "time_ancillary": time,
                    "freq_ancillary": freq,
                    "meta": meta,
                }
            )
        else:
            raise ValueError(f"File format '{self.file_format}' not supported.")

    @staticmethod
    def _read_acq(file_name):
        Q, px, anc = read_acq.decode_file(file_name, progress=False, meta=True)

        freq_anc = {"frequencies": anc.frequencies}
        time_anc = anc.data
        spectra = {"Q": Q, "p0": px[0], "p1": px[1], "p2": px[2]}

        meta = anc.meta
        return spectra, freq_anc, time_anc, meta


@attr.s
class Spectrum(_SpectrumOrResistance):
    """
    Class representing an observed spectrum.

    Standard initialization takes a filename which will be read directly (as long as it
    is in one of the supported formats). Initialization via :func:`from_load` will
    attempt to find a file with the default naming scheme of the database.

    Supported formats: h5, acq, mat

    Examples
    --------
    >>> spec = Spectrum.from_load("Ambient", ".")
    >>> spec.file_format
    h5
    >>> spectra = spec.read()
    """

    supported_formats = ["h5", "acq", "mat"]

    @cached_property
    def _raw_spec(self):
        return FieldSpectrum(self.path)

    @cached_property
    def file_format(self) -> str:
        """The file format of the data to be read."""
        return self._raw_spec.file_format

    @cached_property
    def data(self) -> HDF5RawSpectrum | list[HDF5RawSpectrum]:
        """A view of the data in the file as a HDF5Object.

        If the file is an ACQ file, it will be read completely into memory and cast
        into the same format as a :class:`~h5.HDF5Object` so that the API is the same.

        If the number of files is more than one, `data` will be a list of objects.
        """
        return self._raw_spec.data


@attr.s
class Resistance(_SpectrumOrResistance):
    """An object representing a resistance measurement (and its structure)."""

    supported_formats = ("csv",)

    known_patterns = _SpectrumOrResistance.known_patterns + (
        r"^(?P<load_name>%s)" % _SpectrumOrResistance._loadname_pattern
        + r".(?P<file_format>\w{2,3})$",
    )

    @classmethod
    def from_load(cls, *args, **kwargs):
        classes = super().from_load(*args, **kwargs)
        return classes[0]

    @cached_property
    def file_format(self):
        """The file format of the data to be read."""
        return "csv"

    @classmethod
    def read_csv(cls, path: Path) -> tuple[np.ndarray, dict]:
        with open(path, errors="ignore") as fl:
            if fl.readline().startswith("FLUKE"):
                return cls.read_old_style_csv(path)
            else:
                return cls.read_new_style_csv(path)

    def read(self):
        return self.read_csv(self.path)

    @classmethod
    def read_new_style_csv(cls, path: [str, Path]) -> tuple[np.ndarray, dict]:
        data = np.genfromtxt(
            path,
            skip_header=1,
            delimiter=",",
            dtype=np.dtype(
                [
                    ("date", "S10"),
                    ("time", "S8"),
                    ("lna_voltage", float),
                    ("lna_resistance", float),
                    ("lna_temp", float),
                    ("sp4t_voltage", float),
                    ("sp4t_resistance", float),
                    ("sp4t_temp", float),
                    ("load_voltage", float),
                    ("load_resistance", float),
                    ("load_temp", float),
                    ("room_temp", float),
                ]
            ),
        )
        return data, {}

    @classmethod
    def read_old_style_csv_header(cls, path: Path):
        with open(path, errors="ignore") as fl:
            if not fl.readline().startswith("FLUKE"):
                return {}, 0

            done = False
            out = {}
            nheader_lines = 0
            while not done:
                line = fl.readline()

                if line.startswith("Start Time,") or line.startswith("Max Time,"):
                    names = line.split(",")

                    next_line = fl.readline()
                    nheader_lines += 1
                    values = next_line.split(",")

                    out.update({name: value for name, value in zip(names, values)})

                if line.startswith("1,") or line == "":
                    done = True

                nheader_lines += 1

        return out, nheader_lines

    @classmethod
    def read_old_style_csv(cls, path) -> tuple[np.ndarray, dict]:
        # Weirdly, some old-style files use KOhm, and some just use Ohm.

        # These files have bad encoding, which we can ignore. This means we have to
        # read in the whole thing as text first (while ignoring errors) and construct
        # a StringIO object to pass to genfromtxt.
        header, nheader_lines = cls.read_old_style_csv_header(path)
        nlines = int(header["Total readings"])

        with open(path, errors="ignore") as fl:
            # Get past the header.
            for i in range(nheader_lines):
                next(fl)

            s = StringIO("".join([next(fl) for i in range(nlines - 1)]))

            # Determine whether the file is in KOhm

            def float_from_kohm(x):
                x = x.decode("utf-8")
                kohm = "KOhm" in x
                y = float(x.split(" ")[0])
                return y * 1000 if kohm else y

            data = np.genfromtxt(
                s,
                delimiter=",",
                dtype=np.dtype(
                    [
                        ("reading_num", int),
                        ("sample_resistance", float),
                        ("start_time", "S20"),
                        ("duration", "S9"),
                        ("max_time", "S20"),
                        ("max_resistance", float),
                        ("load_resistance", float),
                        ("min_time", "S20"),
                        ("min_resistance", float),
                        ("description", "S20"),
                        ("end_time", "S22"),
                    ]
                ),
                converters={
                    1: float_from_kohm,
                    5: float_from_kohm,
                    6: float_from_kohm,
                    8: float_from_kohm,
                    10: float_from_kohm,
                },
            )
        return data, {}

    @cached_property
    def _res_and_anc(self):
        return self.read()

    @property
    def resistance(self):
        """The resistance measurement in the file."""
        return self._res_and_anc[0]

    @property
    def ancillary(self):
        """The full raw data from the CSV file."""
        return self._res_and_anc[1]

    @classmethod
    def _get_filename_params_from_contents(cls, path: Path) -> dict:
        meta, _ = cls.read_old_style_csv_header(path)

        if not meta:
            return {}

        start_time = datetime.strptime(meta["Start Time"], "%m/%d/%Y %I:%M:%S %p")

        jd = utils.ymd_to_jd(start_time.year, start_time.month, start_time.day)

        return {
            "hour": start_time.hour,
            "minute": start_time.minute,
            "second": start_time.second,
            "jd": jd,
            "year": start_time.year,
        }


@attr.s
class _SpectraOrResistanceFolder(_DataContainer):
    _run_num: int | dict[str, int] | None = attr.ib(default=None, kw_only=True)
    filetype: str | None = attr.ib(default=None, kw_only=True)

    @cached_property
    def _run_nums(self) -> dict[str, int | None]:
        if type(self._run_num) is int or self._run_num is None:
            return {load: self._run_num for load in LOAD_ALIASES.values()}
        else:
            return self._run_num

    @cached_property
    def _loads(self) -> dict[str, Spectrum | Resistance]:
        loads = {}
        for name, load in LOAD_ALIASES.items():
            try:
                loads[name] = self._content_type.from_load(
                    load, self.path, self._run_nums.get(load, None), self.filetype
                )
            except utils.LoadExistError:
                pass

        return loads

    def __getattr__(self, item):
        if item in self._loads:
            return self._loads[item]

        raise AttributeError(f"{item} does not exist!")

    @cached_property
    def simulators(self) -> dict[str, Spectrum | Resistance]:
        sims = {}
        for name in self.get_simulator_names(self.path):
            sims[name] = self._content_type.from_load(
                name, self.path, self._run_nums.get(name, None), self.filetype
            )
        return sims

    @property
    def load_names(self) -> tuple[str]:
        return tuple(LOAD_ALIASES.keys())

    @property
    def available_load_names(self) -> tuple[str]:
        return tuple(name for name in self.load_names if hasattr(self, name))

    @property
    def run_num(self) -> dict[str, int]:
        """Dictionary of run numbers for each load"""
        try:
            return {k: getattr(self, k)[0].run_num for k in self.available_load_names}
        except TypeError:
            return {k: getattr(self, k).run_num for k in self.available_load_names}

    @classmethod
    def _check_all_files_there(cls, path: Path) -> bool:
        # Just need to check for the loads.
        ok = True
        for name, load in LOAD_ALIASES.items():
            if not path.glob(load + "_*"):
                logger.error(
                    f"{cls.__name__} does not contain any files for load {load}"
                )
                ok = False
        return ok

    @classmethod
    def get_all_load_names(cls, path) -> set[str]:
        """Get all load names found in the Spectra directory"""
        fls = utils.get_active_files(path)
        return {fl.name.split("_")[0] for fl in fls}

    @classmethod
    def get_simulator_names(cls, path) -> set[str]:
        load_names = cls.get_all_load_names(path)
        return {name for name in load_names if name in ANTENNA_SIMULATORS}

    @classmethod
    def _check_file_consistency(cls, path: Path) -> bool:
        fls = utils.get_active_files(path)
        ok = True

        groups = [
            re.search(cls._content_type.pattern, fl.name).groupdict() for fl in fls
        ]

        # Ensure all years are the same
        for fl, group in zip(fls, groups):
            if group["year"] != groups[0]["year"]:
                logger.error(
                    f"All years must be the same in a Spectra folder, but {fl} was not"
                )
                ok = False

        # Ensure days are close-ish
        days = [int(group["day"]) for group in groups]
        if max(days) - min(days) > 30:
            logger.error(f"Observation days are suspiciously far apart for {path}")
            ok = False

        return ok

    def read_all(self):
        """Read all spectra"""
        out = {}
        meta = {}
        for name in self.available_load_names:
            out[name], meta[name] = getattr(self, name).read()
        return out


@attr.s
class Spectra(_SpectraOrResistanceFolder):
    pattern = "Spectra"
    known_patterns = ("spectra",)
    _content_type = Spectrum
    write_pattern = "Spectra"


@attr.s
class Resistances(_SpectraOrResistanceFolder):
    pattern = "Resistance"
    known_patterns = ("resistance",)
    _content_type = Resistance
    write_pattern = "Resistance"


@attr.s
class S1P(_DataFile):
    POSSIBLE_KINDS = [
        "Match",
        "Short",
        "Open",
        "ExternalMatch",
        "ExternalShort",
        "ExternalOpen",
        "External",
        "ReceiverReading",
        "ExternalLoad",
    ]
    pattern = r"^(?P<kind>%s)(?P<repeat_num>\d{2}).s1p$" % ("|".join(POSSIBLE_KINDS))
    write_pattern = "{kind}{repeat_num:>02}.s1p"
    known_patterns = (
        r"^(?P<kind>%s)(?P<repeat_num>\d{1}).s1p$" % ("|".join(POSSIBLE_KINDS)),
        fr"^(?P<kind>{'|'.join(k.lower() for k in POSSIBLE_KINDS)})(?P<repeat_num>\d{2}).s1p$",
        fr"^(?P<kind>{'|'.join(k.lower() for k in POSSIBLE_KINDS)})(?P<repeat_num>\d{1}).s1p$",
        r"^(?P<kind>%s).s1p$" % ("|".join(POSSIBLE_KINDS)),
        fr"^(?P<kind>{'|'.join(k.lower() for k in POSSIBLE_KINDS)}).s1p$",
    )
    known_substitutions = (("Ext_", "External"), ("Int_", ""))  # "Internal"

    @classmethod
    def typestr(cls, name: str) -> str:
        return cls.__name__ + re.match(cls.pattern, name).groupdict()["kind"]

    @property
    def kind(self):
        """The standard of this S1P measurement."""
        return self._match_dict["kind"]

    @property
    def repeat_num(self):
        """The repeat num of this S1P."""
        return self._match_dict["repeat_num"]

    @cached_property
    def s11(self):
        """The S11 measurement in this S1P file.

        Corresponds to :attr:`freq`.
        """
        return self.read(self.path)[0]

    @cached_property
    def freq(self):
        """The frequencies of the S11 measurement in this S1P file.

        Corresponds to :attr:`s11`.
        """
        return self.read(self.path)[1]

    @classmethod
    def _validate_match(cls, match: dict[str, str], filename: str):
        if int(match["repeat_num"]) < 1:
            logger.error(
                f"The file {filename} has a repeat_num ({match['repeat_num']}) less than one"
            )

    @classmethod
    def _get_filename_parameters(cls, dct: dict):
        # If a lower-case kind is passed, use the upper-case version
        out = {"repeat_num": 1}
        if dct.get("kind", None) in (k.lower() for k in cls.POSSIBLE_KINDS):
            dct["kind"] = cls.POSSIBLE_KINDS[
                [k.lower() for k in cls.POSSIBLE_KINDS].index(dct["kind"])
            ]
        return out

    @classmethod
    def read(cls, path_filename):
        d, flag = cls._get_kind(path_filename)
        f = d[:, 0]

        if flag == "DB":
            r = 10 ** (d[:, 1] / 20) * (
                np.cos((np.pi / 180) * d[:, 2]) + 1j * np.sin((np.pi / 180) * d[:, 2])
            )
        elif flag == "MA":
            r = d[:, 1] * (
                np.cos((np.pi / 180) * d[:, 2]) + 1j * np.sin((np.pi / 180) * d[:, 2])
            )
        elif flag == "RI":
            r = d[:, 1] + 1j * d[:, 2]
        else:
            raise Exception("file had no flags set!")

        return r, f / 1e6

    @staticmethod
    def _get_kind(path_filename):
        # identifying the format

        with open(path_filename) as d:
            comment_rows = 0
            flag = None
            for line in d.readlines():
                # checking settings line
                if line.startswith("#"):
                    if "DB" in line or "dB" in line:
                        flag = "DB"
                    if "MA" in line:
                        flag = "MA"
                    if "RI" in line:
                        flag = "RI"

                    comment_rows += 1
                elif line.startswith("!"):
                    comment_rows += 1
                elif flag is not None:
                    break
                else:
                    warnings.warn(
                        f"Non standard line in S11 file {path_filename}: '{line}'\n...Treating as a comment line."
                    )
                    comment_rows += 1

        if flag is None:
            raise OSError(f"The file {path_filename} has incorrect format.")

        #  loading data
        d = np.genfromtxt(path_filename, skip_header=comment_rows)

        return d, flag


@attr.s
class _S11SubDir(_DataContainer):
    STANDARD_NAMES = S1P.POSSIBLE_KINDS
    _content_type = S1P
    write_pattern = "{load_name}{run_num:0>2}"

    repeat_num: int = attr.ib(kw_only=True, converter=int)

    @repeat_num.default
    def _repnum_default(self):
        return self._get_max_repeat_num()

    @cached_property
    def run_num(self) -> int:
        return int(self._match_dict["run_num"])

    @classmethod
    def typestr(cls, name: str) -> str:
        return cls.__name__ + re.match(cls.pattern, name).groupdict()["load_name"]

    @cached_property
    def children(self) -> dict[str, S1P]:
        """Filenames of S1P measurements used in this observation."""
        return {
            name.lower(): S1P(self.path / f"{name}{self.repeat_num:>02}.s1p")
            for name in self.STANDARD_NAMES
        }

    def __getattr__(self, item):
        try:
            return self.children[item]
        except KeyError:
            raise AttributeError(
                f"{item} is not an attribute of {self.__class__.__name__}"
            )

    @cached_property
    def filenames(self) -> tuple[Path]:
        """Filenames of S1P measurements used in this observation."""
        return tuple(val.path for val in self.children.values())

    @property
    def freq(self):
        """Frequencies measured in child S1P files."""
        return self.children["match"].freq

    @property
    def active_contents(self):
        return utils.get_active_files(self.path)

    @classmethod
    def _check_all_files_there(cls, path: Path) -> bool:
        ok = True
        for name in cls.STANDARD_NAMES:
            if not path.glob(name + "??.s1p"):
                logger.error(f"No {name} standard found in {path}")
                ok = False
        return ok

    def _get_max_repeat_num(self) -> int:
        if self.active_contents:
            return max(
                int(re.match(S1P.pattern, fl.name).group("repeat_num"))
                for fl in self.active_contents
            )
        else:
            return 0

    @property
    def max_repeat_num(self) -> int:
        return self._get_max_repeat_num()

    @classmethod
    def _check_file_consistency(cls, path: Path) -> bool:
        return True

    @classmethod
    def _get_filename_parameters(cls, dct: dict):
        out = {}
        if "run_num" not in dct:
            out["run_num"] = 1
        return out


@attr.s
class LoadS11(_S11SubDir):
    STANDARD_NAMES = ["Open", "Short", "Match", "External"]
    pattern = r"(?P<load_name>%s)(?P<run_num>\d{2})$" % (
        "|".join(LOAD_ALIASES.values())
    )
    known_patterns = (
        f"(?P<load_name>{'|'.join(LOAD_MAPPINGS.keys())})$",
        f"(?P<load_name>{'|'.join(LOAD_ALIASES.values())})$",
        r"(?P<load_name>%s)(?P<run_num>\d{1})$" % ("|".join(LOAD_ALIASES.values())),
    )

    known_substitutions = (
        ("AmbientLoad", "Ambient"),
        ("LongCableShort_", "LongCableShorted_"),
    )

    @cached_property
    def load_name(self) -> str:
        return LOAD_ALIASES.inverse.get(
            self._match_dict["load_name"], self._match_dict["load_name"]
        )

    @classmethod
    def _get_filename_parameters(cls, dct: dict):
        out = super()._get_filename_parameters(dct)
        if dct["load_name"] in LOAD_MAPPINGS:
            dct["load_name"] = LOAD_MAPPINGS[dct["load_name"]]
        return out


@attr.s
class AntSimS11(LoadS11):
    pattern = r"(?P<load_name>%s)(?P<run_num>\d{2})$" % (
        "|".join(ANTENNA_SIMULATORS.keys())
    )
    known_patterns = (
        r"(?P<load_name>%s)$" % ("|".join(ANTSIM_REVERSE.keys())),
        r"(?P<load_name>%s)$" % ("|".join(ANTENNA_SIMULATORS.keys())),
    )

    @classmethod
    def _get_filename_parameters(cls, dct: dict) -> dict:
        out = super()._get_filename_parameters(dct)

        if dct["load_name"] in ANTSIM_REVERSE:
            dct["load_name"] = ANTSIM_REVERSE[dct["load_name"]]
        return out


@attr.s
class SwitchingState(_S11SubDir):
    pattern = r"(?P<load_name>SwitchingState)(?P<run_num>\d{2})$"
    known_patterns = ("(?P<load_name>SwitchingState)",)

    STANDARD_NAMES = [
        "Open",
        "Short",
        "Match",
        "ExternalOpen",
        "ExternalShort",
        "ExternalMatch",
    ]
    known_substitutions = (("InternalSwitch", "SwitchingState"),)


@attr.s
class ReceiverReading(_S11SubDir):
    pattern = r"(?P<load_name>ReceiverReading)(?P<run_num>\d{2})$"
    STANDARD_NAMES = ["Open", "Short", "Match", "ReceiverReading"]
    known_substitutions = (("ReceiverReadings", "ReceiverReading"),)
    known_patterns = ("(?P<load_name>ReceiverReading)",)


@attr.s
class S11Dir(_DataContainer):
    """Class representing the entire S11 subdirectory of an observation

    Parameters
    ----------
    path : str or Path
        Top-level directory of the S11 measurements.
    repeat_num : int or dict, optional
        If int, the repeat num of all the standards to use (typically one or two).
        Otherwise, specified as a dict per-load. By default, use the highest repeat
        number available.
    run_num : int or dict, optional
        If int, the run num to use for all loads. If dict, specify which run num
        to use for each load. By default, use the highest run for each available
        load. **Note:** if using this for calibration, the run number must match
        the run number of the spectra and resistance.

    """

    _content_type = {
        **{load: LoadS11 for load in LOAD_ALIASES.values()},
        **{load: LoadS11 for load in LOAD_MAPPINGS.keys()},
        **{
            "SwitchingState": SwitchingState,
            "ReceiverReading": ReceiverReading,
            "InternalSwitch": SwitchingState,  # To catch the old way so it can be fixed.
        },
        **{key: AntSimS11 for key in ANTENNA_SIMULATORS.keys()},
        **{key: AntSimS11 for key in ANTSIM_REVERSE.keys()},
    }
    pattern = "S11"
    known_patterns = ("s11",)
    write_pattern = "S11"

    _repeat_num: int | dict = attr.ib(default=attr.Factory(dict))
    _run_num: int | dict | None = attr.ib(default=attr.Factory(dict))

    @cached_property
    def _run_nums(self) -> dict[str, int]:
        run_nums = {}
        for name in (
            ["switching_state", "receiver_reading"]
            + list(self.available_load_names)
            + list(self.get_simulator_names(self.path))
        ):
            try:
                if isinstance(self._run_num, int):
                    run_nums[name] = self._run_num
                elif isinstance(self._run_num, dict):
                    run_nums[name] = self._run_num.get(
                        name,
                        self._get_highest_run_num(
                            self.path, utils.snake_to_camel(name)
                        ),
                    )
                else:
                    raise ValueError("run_num must be an int or dict.")
            except FileNotFoundError:
                # that's fine, it's probably switching_state or receiver_Reading
                pass

        return run_nums

    @cached_property
    def _repeat_nums(self) -> dict[str, int]:
        if self._repeat_num is None or isinstance(self._repeat_num, int):
            return {
                **{
                    "switching_state": self._repeat_num,
                    "receiver_reading": self._repeat_num,
                },
                **{name: self._repeat_num for name in LOAD_ALIASES.values()},
            }
        else:
            return dict(self._repeat_num)

    @cached_property
    def switching_state(self):
        if "switching_state" in self._run_nums:
            return SwitchingState(
                self.path / f"SwitchingState{self._run_nums['switching_state']:>02}",
                repeat_num=self._repeat_nums.get("switching_state", attr.NOTHING),
            )
        else:
            raise AttributeError("switching_state does not exist")

    @cached_property
    def receiver_reading(self):
        if "receiver_reading" in self._run_nums:
            return ReceiverReading(
                self.path / f"ReceiverReading{self._run_nums['receiver_reading']:>02}",
                repeat_num=self._repeat_nums.get("receiver_reading", attr.NOTHING),
            )
        else:
            raise AttributeError("receiver_reading does not exist")

    @cached_property
    def _loads(self) -> dict[str, LoadS11]:
        return {
            name: LoadS11(
                self.path / f"{LOAD_ALIASES[name]}{self._run_nums[name]:>02}",
                repeat_num=self._repeat_nums.get(
                    LOAD_ALIASES[name], self._repeat_nums.get(name, attr.NOTHING)
                ),
            )
            for name in self.available_load_names
        }

    @cached_property
    def simulators(self) -> dict[str, AntSimS11]:
        return {
            name: AntSimS11(
                self.path / f"{name}{self._run_nums[name]:>02}",
                repeat_num=self._repeat_nums.get(name, attr.NOTHING),
            )
            for name in self.get_simulator_names(self.path)
        }

    def __getattr__(self, item):
        if item in self._loads:
            return self._loads[item]

        raise AttributeError(f"{item} does not exist!")

    @property
    def available_load_names(self):
        return self.get_available_load_names(self.path)

    @classmethod
    def get_available_load_names(cls, path):
        fls = utils.get_active_files(path)
        return {
            LOAD_ALIASES.inverse[fl.name[:-2]]
            for fl in fls
            if any(fl.name.startswith(k) for k in LOAD_ALIASES.inverse)
        }

    @property
    def load_names(self) -> tuple[str]:
        return tuple(LOAD_ALIASES.keys())

    @property
    def repeat_num(self) -> dict[str, int]:
        """Dictionary specifying run numbers for each load."""
        out = {k: getattr(self, k).repeat_num for k in list(self.available_load_names)}

        try:
            out["switching_state"] = self.switching_state.repeat_num
        except AttributeError:
            pass

        try:
            out["receiver_reading"] = self.receiver_reading.repeat_num
        except AttributeError:
            pass

        out.update({k: v.repeat_num for k, v in self.simulators.items()})
        return out

    @property
    def run_num(self) -> dict[str, int]:
        out = {k: getattr(self, k).run_num for k in list(self.available_load_names)}
        try:
            out["switching_state"] = self.switching_state.run_num
        except AttributeError:
            pass

        try:
            out["receiver_reading"] = self.receiver_reading.run_num
        except AttributeError:
            pass

        out.update({k: v.run_num for k, v in self.simulators.items()})
        return out

    @classmethod
    def _get_highest_run_num(cls, path, kind) -> int:
        fls = utils.get_active_files(path)
        fls = [fl for fl in fls if kind in str(fl)]
        if not fls:
            raise FileNotFoundError(f"No S11 measurements found for {kind}")

        run_nums = [int(str(fl)[-2:]) for fl in fls]
        return max(run_nums)

    def get_highest_run_num(self, kind: str) -> int:
        """Get the highest repeat number for this kind."""
        return self._get_highest_run_num(self.path, kind)

    @classmethod
    def _check_all_files_there(cls, path: Path) -> bool:
        ok = True
        for load in LOAD_ALIASES.values():
            if not path.glob(load):
                logger.error(f"No {load} S11 directory found!")
                ok = False

        for other in ["SwitchingState", "ReceiverReading"]:
            if not path.glob(other + "??"):
                logger.error(f"No {other} S11 directory found!")
                ok = False
        return ok

    @classmethod
    def _check_file_consistency(cls, path: Path) -> bool:
        simulators = cls.get_simulator_names(path)
        if simulators:
            logger.info(
                f"Found the following Antenna Simulators in S11: {','.join(simulators)}"
            )
        else:
            logger.info("No Antenna Simulators in S11.")
        return True

    @classmethod
    def get_simulator_names(cls, path) -> set[str]:
        fls = utils.get_active_files(path)
        return {
            fl.name[:-2]
            for fl in fls
            if any(fl.name.startswith(k) for k in ANTENNA_SIMULATORS)
        }


@attr.s
class CalibrationObservation(_DataContainer):
    """
    A full set of data required to calibrate field observations.

    Incorporates several lower-level objects, such as :class:`Spectrum`,
    :class:`Resistance` and :class:`S1P` in a seamless way.

    Parameters
    ----------
    path : str or Path
        The path (absolute or relative to current directory) to the top level
        directory of the observation. This should look something like
        ``Receiver01_2020_01_01_040_to_200MHz/``.
    run_num : int or dict, optional
        If an integer, the run number to use for all loads. If None, by default
        uses the last run for each load. If a dict, it should specify the
        run number for each load.
    repeat_num : int or dict, optional
        If an integer, the repeat number to use for all S11 standards measurements,
        for all loads. If None, by default uses the last repeat (typically 02) for
        each load. If a dict, it should specify the repeat number for each load.
    include_previous : bool, optional
        Whether to by default include the previous observation in the same directory
        to supplement the current one if parts are missing.
    compile_from_def : bool, optional
        Whether to attempt compiling a virtual observation from a ``definition.yaml``
        inside the observation directory. This is the default behaviour, but can
        be turned off to enforce that the current directory should be used directly.
    """

    pattern = re.compile(
        r"^Receiver(?P<rcv_num>\d{2})_(?P<temp>\d{2})C_(?P<year>\d{4})_(?P<month>\d{2})_("
        r"?P<day>\d{2})_(?P<freq_low>\d{3})_to_(?P<freq_hi>\d{3})MHz$"
    )

    known_patterns = (
        (
            r"^Receiver(\d{1,2})_(?P<temp>\d{2})C_(\d{4})_(\d{1,2})_(\d{1,2})_(\d{2,3})_"
            r"to_(\d{2,3})MHz$"
        ),
        (
            r"Receiver(?P<rcv_num>\d{2})_(?P<temp>\d{2})C_(?P<year>\d{4})_(?P<month>\d{2})_"
            r"(?P<day>\d{2})_(?P<freq_low>\d{3})"
            r"_to_(?P<freq_hi>\d{3})_MHz$"
        ),
    )
    write_pattern = (
        "Receiver{rcv_num:0>2}_{temp:>02}C_{year:>04}_{month:0>2}_{day:0>2}_"
        "{freq_low:0>3}_to_{freq_hi:0>3}MHz"
    )

    _content_type = {
        "S11": S11Dir,
        "Spectra": Spectra,
        "Resistance": Resistances,
        "spectra": Spectra,
        "resistance": Resistances,
        "s11": S11Dir,
    }

    _run_num: int | dict = attr.ib(default=attr.Factory(dict))
    _repeat_num: int | dict = attr.ib(default=attr.Factory(dict))
    spectra_kwargs: dict = attr.ib(default=attr.Factory(dict))
    s11_kwargs: dict = attr.ib(default=attr.Factory(dict))
    resistance_kwargs: dict = attr.ib(default=attr.Factory(dict))
    original_path: Path = attr.ib(converter=Path)
    _tmpdir: Path | None = attr.ib(default=None)

    @cached_property
    def definition(self) -> dict:
        # Read the definition file, and combine other observations into a single
        # temporary directory if they exist (otherwise, just symlink this full directory)
        # Note that we need to keep the actual _tmpdir object around otherwise it gets
        # cleaned up!

        definition = self.check_definition(self.original_path)

        if definition.get("entirely_invalid", False):
            logger.warning(
                f"Observation {self.original_path} is marked as invalid -- "
                f"proceed with caution! Reason: '{self.definition['entirely_invalid']}'"
            )

        return definition

    @original_path.default
    def _original_path_default(self):
        # the original input path, so we have access to it later. Otherwise
        # we might just have a temporary directory.
        return self._path

    @classmethod
    def from_def(
        cls, path: str | Path, include_previous: bool = True, **kwargs
    ) -> CalibrationObservation:
        tmpdir, name = cls.compile_obs_from_def(path, include_previous)
        new_path = tmpdir / name
        return cls(new_path, tmpdir=tmpdir, original_path=path, **kwargs)

    @property
    def receiver_num(self) -> int:
        return int(self._match_dict["rcv_num"])

    @property
    def ambient_temp(self) -> int:
        return int(self._match_dict["temp"])

    @property
    def year(self) -> int:
        return int(self._match_dict["year"])

    @property
    def month(self) -> int:
        return int(self._match_dict["month"])

    @property
    def day(self) -> int:
        return int(self._match_dict["day"])

    @property
    def freq_low(self) -> int:
        return int(self._match_dict["freq_low"])

    @property
    def freq_high(self) -> int:
        return int(self._match_dict["freq_hi"])

    @cached_property
    def spectra(self) -> Spectra:
        return Spectra(
            self.path / "Spectra",
            run_num=self._run_num,
            **self.spectra_kwargs,
        )

    @cached_property
    def resistance(self) -> Resistances:
        return Resistances(
            self.path / "Resistance",
            run_num=self._run_num,
            **self.resistance_kwargs,
        )

    @cached_property
    def s11(self) -> S11Dir:
        return S11Dir(
            self.path / "S11",
            run_num=self._run_num,
            repeat_num=self._repeat_num,
            **self.s11_kwargs,
        )

    @cached_property
    def simulator_names(self):
        return self.get_simulator_names(self.path)

    @classmethod
    def from_observation_yaml(cls, obs_yaml: [str, Path]):
        """Create a CalibrationObservation from a specific YAML format."""
        obs_yaml = Path(obs_yaml)
        assert obs_yaml.exists(), f"{obs_yaml} does not exist!"

        with open(obs_yaml) as fl:
            obs_yaml_data = yaml.load(fl, Loader=yaml.FullLoader)

        root = obs_yaml_data["root"]
        root = obs_yaml.parent.absolute() if not root else Path(root).absolute()
        assert (
            root.exists()
        ), f"The root {root} specified in the observation does not exist."

        files = obs_yaml_data["files"]
        meta = obs_yaml_data["meta"]
        cls._check_yaml_files(files, root)

        tmpdir = tempfile.TemporaryDirectory()

        sympath = Path(tmpdir.name) / cls.write_pattern.format(**meta)
        sympath.mkdir(parents=True)

        # Make top-level directories
        spec = sympath / "Spectra"
        s11 = sympath / "S11"
        res = sympath / "Resistance"
        spec.mkdir()
        s11.mkdir()
        res.mkdir()

        # Link all Spectra and Resistance files.
        for key, thing in zip(["spectra", "resistance"], [spec, res]):
            for kind, kind_files in files[key].items():
                these_files = sum((list(root.glob(fl)) for fl in kind_files), [])
                for fl in these_files:
                    (thing / fl.name).symlink_to(root / fl)

        # Symlink the S11 files.
        s11_run_nums = {}
        rep_nums = {}
        for key, (direc, run_num) in files["s11"].items():
            direc = Path(root / direc)
            syms11 = s11 / direc.name
            syms11.mkdir()

            if key == "receiver":
                rep_nums["receiver_reading"] = int(str(direc)[-2:])
            elif key == "switch":
                rep_nums["switching_state"] = int(str(direc)[-2:])

            these_files = direc.glob(f"*{run_num:>02}.?1?")
            for fl in these_files:
                (syms11 / fl.name).symlink_to(direc / fl.name)

            s11_run_nums[key] = run_num

        return cls(
            sympath,
            run_num={"S11": s11_run_nums},
            repeat_num=rep_nums,
            fix=False,
            tmpdir=tmpdir,
        )

    @classmethod
    def _check_yaml_files(cls, files: dict, root: Path):
        """Check goodness of 'files' key in an observation yaml."""
        for key in ["spectra", "resistance", "s11"]:
            assert key in files, f"{key} must be in observation YAML 'files'"
            for key2 in ["open", "short", "hot_load", "ambient"]:
                assert (
                    key2 in files[key]
                ), f"{key2} must be in observation YAML 'files.{key}'"

                if key != "s11":
                    for fl in files[key][key2]:
                        assert (
                            len(list(root.glob(fl))) > 0
                        ), f"File '{root / fl}' included at files.{key}.{key2} does not exist or match any glob patterns."
                else:
                    fl = files[key][key2][0]
                    assert (
                        root / fl
                    ).exists(), f"Directory '{root / fl}' included at files.{key}.{key2} does not exist."

            if key == "s11":
                for key2 in ["receiver", "switch"]:
                    assert (
                        key2 in files[key]
                    ), f"{key2} must be in observation YAML 'files.{key}'. Available: {list(files[key].keys())}"

                    assert (
                        root / files[key][key2][0]
                    ).exists(), f"Directory '{root /files[key][key2][0]}' included at files.{key}.{key2} does not exist."

    @classmethod
    def check_definition(cls, path: Path) -> dict:
        """Check the associated definition.yaml file within an observation."""

        definition_file = path / "definition.yaml"

        # Read in the definition file (if it exists)
        if not definition_file.exists():
            return {}

        with open(definition_file) as fl:
            definition = yaml.load(fl, Loader=yaml.FullLoader) or {}

        allowed_keys = {
            "root_obs_dir": str,
            "entirely_invalid": str,
            "include": list,
            "prefer": list,
            "invalid": list,
            "measurements": {
                "resistance_m": dict(
                    tuple((i, float) for i in range(99))
                    + tuple((f"{i:02}", float) for i in range(99))
                ),
                "resistance_f": dict(
                    tuple((i, float) for i in range(99))
                    + tuple((f"{i:02}", float) for i in range(99))
                ),
            },
            "defaults": {"run": dict, "repeat": dict},
            "purpose": str,
            "history": list,
        }

        def _check_grp(defn, allowed):
            for k, v in defn.items():
                if k not in allowed:
                    logger.warning(
                        f"Key {k} found in definition.yaml, but is not a known keyword."
                    )
                elif isinstance(allowed[k], dict):
                    # Recurse into sub-dictionaries.
                    _check_grp(v, allowed[k])
                elif not isinstance(v, allowed[k]):
                    logger.error(
                        f"Key {k} has wrong type in definition.yaml. Should be {allowed[k]}, got {type(v)}."
                    )

        _check_grp(definition, allowed_keys)
        return definition

    @classmethod
    def _check_self(cls, path: Path, **kwargs):
        path = path.absolute()

        # Warn if this is an invalid observation entirely. Also, we don't check the
        # observation then, as it's annoyingly difficult.
        if path.parent.suffix in [".invalid", ".old"]:
            logger.warning(
                f"Observation {path.parent.name} is marked as {path.parent.suffix} -- "
                f"proceed with caution!"
            )
            return path, None

        return super()._check_self(path, **kwargs)

    @classmethod
    def _validate_match(cls, match: dict[str, str], filename: str):
        groups = match
        if int(groups["rcv_num"]) < 1:
            logger.error(f"Unknown receiver number: {groups['rcv_num']}")
        if not (2010 <= int(groups["year"]) <= 2030):
            logger.error(f"Unknown year: {groups['year']}")
        if not (1 <= int(groups["month"]) <= 12):
            logger.error(f"Unknown month: {groups['month']}")
        if not (1 <= int(groups["day"]) <= 31):
            logger.error(f"Unknown day: {groups['day']}")
        if not (1 <= int(groups["freq_low"]) <= 300):
            logger.error(f"Low frequency is weird: {groups['freq_low']}")
        if not (1 <= int(groups["freq_hi"]) <= 300):
            logger.error(f"High frequency is weird: {groups['freq_hi']}")
        if not int(groups["freq_low"]) < int(groups["freq_hi"]):
            logger.error(
                f"Low frequency > High Frequency: {groups['freq_low']} > {groups['freq_hi']}"
            )
        if not (0 < int(groups["temp"]) < 100):
            logger.error(
                f"Ambient temperature out of range (0 - 100): {groups['temp']}"
            )

        logger.info("Calibration Observation Metadata:")
        for k, v in groups.items():
            logger.info(f"\t{k}: {v}")

    @classmethod
    def path_to_datetime(cls, path: [str, Path]):
        pre_level = logger.getEffectiveLevel()
        logger.setLevel(39)
        try:
            path, match = cls.check_self(path, fix=False)
        except Exception as e:
            raise e
        finally:
            logger.setLevel(pre_level)

        if match:
            return datetime(int(match["year"]), int(match["month"]), int(match["day"]))
        else:
            raise utils.FileStructureError("The path is not valid for an Observation.")

    @classmethod
    def _check_all_files_there(cls, path: Path) -> bool:
        ok = True
        for folder in ["S11", "Spectra", "Resistance"]:
            if not (path / folder).exists():
                logger.warning(f"No {folder} folder in observation!")
                ok = False
        return ok

    @classmethod
    def _check_file_consistency(cls, path: Path) -> bool:
        # checks whether simulators are the same between each.
        cls.get_simulator_names(path)
        return True

    @classmethod
    def get_simulator_names(cls, path: [str, Path]):
        # Go through the subdirectories and check their simulators
        path = Path(path)
        dct = {
            name: tuple(
                sorted(cls._content_type[name].get_simulator_names(path / name))
            )
            for name in ["Spectra", "S11", "Resistance"]
            if (path / name).exists()
        }

        # If any list of simulators is not the same as the others, make an error.
        if len(set(dct.values())) != 1:
            logger.warning(
                f"Antenna Simulators do not match in all subdirectories. Got {dct}"
            )
            names = [
                name
                for name in ANTENNA_SIMULATORS
                if all(name in val for val in dct.values())
            ]
        else:
            names = list(dct.values())[0]

        return set(names)

    def read_all(self):
        """Read all spectra and resistance files into memory. Usually a bad idea."""
        self.spectra.read_all()
        self.resistance.read_all()

    @classmethod
    def get_base_files(cls, path: Path, with_notes=False) -> list[Path]:
        """Get a list of valid files in this observation.

        Takes into account the definition.yaml if it exists.
        """
        definition = cls.check_definition(path)

        invalid = []
        for pattern in definition.get("invalid", []):
            invalid.extend(path.glob(pattern))

        other_ignores = ["definition.yaml"]
        if not with_notes:
            other_ignores.append("Notes.txt")

        # We'll get everything in this subtree, except those marked invalid.
        return utils.get_file_list(
            path,
            filter=lambda x: x.suffix not in [".invalid", ".old"]
            and x.name not in other_ignores
            and x.parent.name != "outputs",
            ignore=invalid,
        )

    @classmethod
    def compile_obs_from_def(
        cls, path: Path, include_previous=True
    ) -> tuple[Path, str]:
        """Make a tempdir containing pointers to relevant files built from a definition.

        Takes a definition file (YAML format) from a particular Observation, and uses
        the ``include`` and ``prefer`` sections to generate a full list of files from any
        number of Observations that make up a single full observation. Will only include
        a single file of each kind.

        Parameters
        ----------
        path : Path
            The path (absolute or relative to current directory) to the observation (not
            the definition file).
        include_previous : bool, optional
            Whether to by default "include" the previous observation (if any can be
            found). This means that observation will be used to supplement this one if
            this is incomplete.

        Returns
        -------
        TemporaryDirectory :
            A temp directory in which there will be a directory of the same name as
            ``path``, and under which will be a view into a "full" observation, compiled
            from the definition. Each file in the directory will be a symlink.
            Note that this return variable must be kept alive or the directory will be
            cleaned up.
        name : str
            The name of the observation (i.e. the directory inside the temporary direc).
        """
        path = Path(path).absolute()
        obs_name = path.name

        assert path.exists(), f"{path} does not exist!"
        definition = cls.check_definition(path)

        # Now include files from other observations if they don't already exist.
        root_obs = definition.get("root_obs_dir", None)

        if root_obs is None:
            root_obs = path.parent
        else:
            # Root observation directory should be relative to the definition file.
            if not Path(root_obs).is_absolute():
                root_obs = (path / root_obs).resolve()

        files = {fl.relative_to(path.parent): fl for fl in cls.get_base_files(path)}

        file_parts = {
            fl.relative_to(obs_name): cls.match_path(fl, root=path.parent)
            for fl in files
        }
        # Actually need files to *not* have the top-level name
        files = {fl.relative_to(obs_name): fl_abs for fl, fl_abs in files.items()}

        def _include_extra(roots, prefer):
            for inc_path in roots:
                # Need to get this root_obs if inc_path is absolute, because we need
                # to know where the observation starts (in the path)
                if inc_path.is_absolute():
                    for indx, part in enumerate(inc_path.parts[::-1]):
                        if cls.pattern.search(part):
                            break
                    else:
                        raise ValueError(
                            f"Can't find an observation root in {inc_path}"
                        )

                    this_root_obs = inc_path.parents[indx]
                    inc_path = inc_path.relative_to(this_root_obs)
                else:
                    this_root_obs = root_obs
                    inc_path = this_root_obs / inc_path
                this_obs_name = inc_path.relative_to(this_root_obs).parts[0]

                # Get all non-invalid files in the other observation.
                inc_files = cls.get_base_files(inc_path)

                # Get the defining classes for each file
                inc_file_parts = {
                    fl: cls.match_path(
                        fl.relative_to(this_root_obs), root=this_root_obs
                    )
                    for fl in inc_files
                }

                new_file_parts = {}
                # Check if the defining classes are the same as any already in there.
                for inc_fl, kinds in inc_file_parts.items():
                    if prefer or not any(kinds == k for k in file_parts.values()):
                        if prefer:
                            # First delete the thing that's already there
                            for k, v in list(file_parts.items()):
                                if v == kinds:
                                    del file_parts[k]
                                    del files[k]

                        files[
                            inc_fl.relative_to(this_root_obs / this_obs_name)
                        ] = inc_fl
                        new_file_parts[inc_fl.relative_to(this_root_obs)] = kinds

                # Updating the file parts after the full loop means that we can add
                # multiple files of the same kind (eg. with different run_num) from a
                # single observation, but only if they didn't exist in a previous obs.
                file_parts.update(new_file_parts)

        default_includes = []
        if include_previous:
            # Look for a previous definition in the root observation directory.
            potential_obs = root_obs.glob(obs_name.split("_")[0] + "_*")
            potential_obs = sorted(
                str(p.name) for p in list(potential_obs) + [path.parent]
            )
            if len(potential_obs) > 1:
                indx = potential_obs.index(obs_name) - 1
                if indx >= 0:
                    default_includes.append(potential_obs[indx])

        include = [Path(x) for x in definition.get("include", default_includes)]
        prefer = [Path(x) for x in definition.get("prefer", [])]
        _include_extra(include, prefer=False)
        _include_extra(prefer, prefer=True)

        stuff = f"{path}_{include_previous}"
        if (path / "definition.yaml").exists():
            # Now make a full symlink directory with these files.
            with open(path / "definition.yaml") as fl:
                stuff += fl.read()
        hsh = hash(stuff)
        dirname = f"calobs_{hsh}"

        symdir = Path(tempfile.gettempdir()) / dirname

        if not symdir.exists():
            symdir.mkdir()
            utils.make_symlink_tree(files, symdir, obs_name)

        return symdir, obs_name

    @classmethod
    def match_path(
        cls, path: Path, root: Path = Path()
    ) -> tuple[_DataFile | _DataContainer]:
        """Give a path relative to the root, determine its describing class.

        Examples
        --------
        >>> CalibrationObservation.match_path('Spectra')
        >>> (Spectra, )
        """
        structure = {
            CalibrationObservation: {
                Spectra: (Spectrum,),
                Resistances: (Resistance,),
                S11Dir: {
                    LoadS11: (S1P,),
                    AntSimS11: (S1P,),
                    SwitchingState: (S1P,),
                    ReceiverReading: (S1P,),
                },
            }
        }

        pre_level = logger.level
        if logger.handlers:
            logger.handlers[0].setLevel(100)  # Temporarily disable stdout handler

        # Add a string buffer handler that can capture the error messages.
        msg_buffer = StringIO()
        handler = logging.StreamHandler(msg_buffer)
        logger.addHandler(handler)

        _strc = structure

        # Get parts of the path, but keep the top-level and the '25C' together
        path_parts = list(path.parts)

        try:
            parts = ()
            full_part = root
            for part in path_parts:
                full_part = Path(full_part) / Path(part)

                for thing in _strc:
                    pth, match = thing.check_self(part, fix=False)

                    if match is not None:
                        parts = parts + (thing.typestr(part),)

                        if isinstance(_strc, dict):
                            _strc = _strc[thing]

                        # Rewind buffer to start afresh with captured errors.
                        msg_buffer.truncate(0)
                        msg_buffer.seek(0)
                        break
                else:
                    raise ValueError(
                        f"path {path} does not seem to point to a known kind of object. "
                        f"Stuck on {part}. Errors/comments received:\n\n{msg_buffer.getvalue()}"
                    )
        except ValueError as e:
            raise e
        finally:
            logger.removeHandler(handler)
            if logger.handlers:
                logger.handlers[0].setLevel(pre_level)

        return parts

    @property
    def run_num(self) -> dict[str, int]:
        """Dictionary specifying run numbers for each component"""
        return self.s11.run_num

    @property
    def list_of_files(self):
        """A list of all data files used in this observation."""
        fls = []
        for name in self.s11.load_names:
            fls += list(getattr(self.s11, name).filenames)

        for name in self.spectra.load_names:
            fls += [x.path for x in getattr(self.spectra, name)]
            fls.append(getattr(self.resistance, name).path)

        return sorted(fl.resolve() for fl in fls)
