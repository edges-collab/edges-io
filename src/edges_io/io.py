"""
This module defines the overall file structure and internal contents of the
calibration observations. It does *not* implement any algorithms/methods on that data,
making it easier to separate the algorithms from the data checking/reading.
"""
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


class _SpectrumOrResistance(_DataFile):
    load_pattern = "|".join(LOAD_ALIASES.values())
    antsim_pattern = "|".join(ANTENNA_SIMULATORS.keys())
    _antsim_rev_pattern = "|".join(ANTSIM_REVERSE.keys())
    _load_rev_pattern = "|".join(LOAD_MAPPINGS.keys())
    _loadname_pattern = (
        f"{load_pattern}|{antsim_pattern}|{_antsim_rev_pattern}|{_load_rev_pattern}"
    )

    pattern = (
        r"(?P<load_name>%s|%s)" % (load_pattern, antsim_pattern)
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
            "^(?P<load_name>{})".format(_loadname_pattern)
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
    def _validate_match(cls, match: Dict[str, str], filename: str):
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
        direc: [str, Path],
        run_num: Optional[int] = None,
        filetype: Optional[str] = None,
    ):
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

        files = list(direc.glob(f"{load}_??_????_???_??_??_??_lab.*"))

        if not files:
            raise ValueError(
                f"No Spectrum files found for load {load}. Available spectrum files: {list(direc.glob('*'))}"
            )

        if filetype:
            files = [fl for fl in files if fl.endswith("." + filetype)]
            if not files:
                raise ValueError(
                    f"No files exist for the load {load} with filetype {filetype} on the path: {direc}"
                )

        else:
            # Use any format so long as it is supported
            restricted_files = []
            for ftype in cls.supported_formats:
                restricted_files = [fl for fl in files if fl.suffix == ("." + ftype)]
                if restricted_files:
                    break

            if not restricted_files:
                raise ValueError(
                    f"No files exist for the load {load} for any filetype on that path: {direc}."
                    f"Found files: {list(files)}."
                )

            files = restricted_files

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
    def days(self) -> int:
        return int(self._match_dict["day"])

    @cached_property
    def load_name(self):
        return LOAD_ALIASES.inverse.get(
            self._match_dict["load_name"], self._match_dict["load_name"]
        )

    @cached_property
    def hours(self):
        """List of integer hours (one per file) at which data acquisition was begun"""
        return int(self._match_dict["hour"])

    @cached_property
    def minutes(self):
        """List of integer minutes (one per file) at which data acquisition was begun"""
        return int(self._match_dict["minute"])

    @cached_property
    def seconds(self):
        """List of integer seconds (one per file) at which data acquisition was begun"""
        return int(self._match_dict["second"])

    def __eq__(self, other):
        return (
            other.__class__.__name__ == self.__class__.__name__
            and self.load_name == other.load_name
        )


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
    def file_format(self) -> str:
        """The file format of the data to be read."""
        return self.path.suffix[1:]

    @cached_property
    def data(self) -> [HDF5RawSpectrum, List[HDF5RawSpectrum]]:
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


class Resistance(_SpectrumOrResistance):
    """An object representing a resistance measurement (and its structure)."""

    supported_formats = ("csv",)

    known_patterns = _SpectrumOrResistance.known_patterns + (
        r"^(?P<load_name>%s)" % _SpectrumOrResistance._loadname_pattern
        + r".(?P<file_format>\w{2,3})$",
    )

    def __init__(self, *args, store_data=True, **kwargs):
        super().__init__(*args, **kwargs)
        self.store_data = store_data

    @classmethod
    def from_load(cls, *args, **kwargs):
        classes = super().from_load(*args, **kwargs)
        return classes[0]

    @cached_property
    def file_format(self):
        """The file format of the data to be read."""
        return "csv"

    @classmethod
    def read_csv(cls, path: Path) -> Tuple[np.ndarray, Dict]:
        with open(path, "r", errors="ignore") as fl:
            if fl.readline().startswith("FLUKE"):
                return cls.read_old_style_csv(path)
            else:
                return cls.read_new_style_csv(path)

    def read(self):
        try:
            return self._data, self._meta
        except AttributeError:
            data, meta = self.read_csv(self.path)

            if self.store_data:
                self._data = data
                self._meta = meta

            return data, meta

    @classmethod
    def read_new_style_csv(cls, path: [str, Path]) -> Tuple[np.ndarray, Dict]:
        data = np.genfromtxt(
            path,
            skip_header=1,
            delimiter=",",
            dtype=np.dtype(
                [
                    ("date", "S10"),
                    ("time", "S8"),
                    ("lna_voltage", np.float),
                    ("lna_resistance", np.float),
                    ("lna_temp", np.float),
                    ("sp4t_voltage", np.float),
                    ("sp4t_resistance", np.float),
                    ("sp4t_temp", np.float),
                    ("load_voltage", np.float),
                    ("load_resistance", np.float),
                    ("load_temp", np.float),
                    ("room_temp", np.float),
                ]
            ),
        )
        return data, {}

    @classmethod
    def read_old_style_csv_header(cls, path: Path):
        with open(path, "r", errors="ignore") as fl:
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
    def read_old_style_csv(cls, path) -> Tuple[np.ndarray, Dict]:
        # Weirdly, some old-style files use KOhm, and some just use Ohm.

        # These files have bad encoding, which we can ignore. This means we have to
        # read in the whole thing as text first (while ignoring errors) and construct
        # a StringIO object to pass to genfromtxt.
        header, nheader_lines = cls.read_old_style_csv_header(path)
        nlines = int(header["Total readings"])

        with open(path, "r", errors="ignore") as fl:
            # Get past the header.
            for i in range(nheader_lines):
                next(fl)

            s = StringIO("".join([next(fl) for i in range(nlines - 1)]))

            # Determine whether the file is in KOhm
            kohm = "KOhm" in s.readline()
            s.seek(0)

            def float_from_kohm(x):
                y = float(x.decode("utf-8").split(" ")[0])
                return y * 1000 if kohm else y

            data = np.genfromtxt(
                s,
                delimiter=",",
                dtype=np.dtype(
                    [
                        ("reading_num", np.int),
                        ("sample_resistance", np.float),
                        ("start_time", "S20"),
                        ("duration", "S9"),
                        ("max_time", "S20"),
                        ("max_resistance", np.float),
                        ("load_resistance", np.float),
                        ("min_time", "S20"),
                        ("min_resistance", np.float),
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

    @property
    def resistance(self):
        """The resistance measurement in the file.

        Note that this is only cached in memory if `store_data` is True, otherwise
        the data is re-read from disk each time `resistance` is accessed.
        """
        if not self.store_data:
            warnings.warn(
                "'resistance' is not being cached -- using it reads the resistance file each time it is accessed!"
            )

        return self.read()[0]

    @property
    def ancillary(self):
        """The full raw data from the CSV file.

        Note that this is only cached in memory if `store_data` is True, otherwise
        the data is re-read from disk each time `resistance` is accessed.
        """
        if not self.store_data:
            warnings.warn(
                "'ancillary' is not being cached -- using it reads the resistance file each time it is accessed!"
            )

        return self.read()[1]

    @classmethod
    def _get_filename_params_from_contents(cls, path: Path) -> Dict:
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


class _SpectraOrResistanceFolder(_DataContainer):
    def __init__(
        self,
        path: [str, Path],
        *,
        run_num: Optional[Union[int, Dict[str, int]]] = None,
        filetype: Optional[str] = None,
        **kwargs,
    ):
        """Collection of spectra in an observation"""
        super().__init__(path, **kwargs)

        if type(run_num) is int or run_num is None:
            run_nums = {load: run_num for load in LOAD_ALIASES.values()}
        else:
            run_nums = run_num

        for name, load in LOAD_ALIASES.items():
            setattr(
                self,
                name,
                self._content_type.from_load(
                    load, self.path, run_nums.get(load, None), filetype
                ),
            )

        # Populate simulators.
        self.simulators = {}
        for name in self.get_simulator_names(self.path):
            self.simulators[name] = self._content_type.from_load(
                name, self.path, run_nums.get(name, None), filetype
            )

    @property
    def load_names(self):
        return tuple(LOAD_ALIASES.keys())

    @property
    def run_num(self):
        """Dictionary of run numbers for each load"""
        try:
            return {k: getattr(self, k)[0].run_num for k in LOAD_ALIASES}
        except TypeError:
            return {k: getattr(self, k).run_num for k in LOAD_ALIASES}

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
    def get_all_load_names(cls, path):
        """Get all load names found in the Spectra directory"""
        fls = utils.get_active_files(path)
        return {fl.name.split("_")[0] for fl in fls}

    @classmethod
    def get_simulator_names(cls, path):
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
        for name in LOAD_ALIASES:
            out[name], meta[name] = getattr(self, name).read()
        return out

    def __eq__(self, other):
        return self.__class__.__name__ == other.__class__.__name__


class Spectra(_SpectraOrResistanceFolder):
    pattern = "Spectra"
    known_patterns = ("spectra",)
    _content_type = Spectrum
    write_pattern = "Spectra"


class Resistances(_SpectraOrResistanceFolder):
    pattern = "Resistance"
    known_patterns = ("resistance",)
    _content_type = Resistance
    write_pattern = "Resistance"


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
    pattern = r"^(?P<kind>%s)(?P<run_num>\d{2}).s1p$" % ("|".join(POSSIBLE_KINDS))
    write_pattern = "{kind}{run_num:>02}.s1p"
    known_patterns = (
        r"^(?P<kind>%s)(?P<run_num>\d{1}).s1p$" % ("|".join(POSSIBLE_KINDS)),
        fr"^(?P<kind>{'|'.join(k.lower() for k in POSSIBLE_KINDS)})(?P<run_num>\d{2}).s1p$",
        fr"^(?P<kind>{'|'.join(k.lower() for k in POSSIBLE_KINDS)})(?P<run_num>\d{1}).s1p$",
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
    def run_num(self):
        """The run num of this S1P."""
        return self._match_dict["run_num"]

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
    def _validate_match(cls, match: Dict[str, str], filename: str):
        if int(match["run_num"]) < 1:
            logger.error(
                f"The file {filename} has a run_num ({match['run_num']}) less than one"
            )

    @classmethod
    def _get_filename_parameters(cls, dct: dict):
        # If a lower-case kind is passed, use the upper-case version
        out = {"run_num": 1}
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

        with open(path_filename, "r") as d:
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
            raise IOError(f"The file {path_filename} has incorrect format.")

        #  loading data
        d = np.genfromtxt(path_filename, skip_header=comment_rows)

        return d, flag

    def __eq__(self, other):
        return (
            self.__class__.__name__ == other.__class__.__name__
            and self.kind == other.kind
        )


class _S11SubDir(_DataContainer):
    STANDARD_NAMES = S1P.POSSIBLE_KINDS
    _content_type = S1P
    write_pattern = "{load_name}{repeat_num:0>2}"

    @classmethod
    def typestr(cls, name: str) -> str:
        return cls.__name__ + re.match(cls.pattern, name).groupdict()["load_name"]

    def __init__(self, path, *, run_num=None, **kwargs):
        super().__init__(path, **kwargs)

        self.run_num = run_num or self._get_max_run_num()
        self.repeat_num = int(self._match_dict["repeat_num"])

    @cached_property
    def children(self) -> Dict[str, S1P]:
        """Filenames of S1P measurements used in this observation."""
        return {
            name.lower(): S1P(self.path / f"{name}{self.run_num:>02}.s1p")
            for name in self.STANDARD_NAMES
        }

    def __getattr__(self, item):
        try:
            return super().__getattr__(item)
        except AttributeError as e:
            try:
                return self.children[item]
            except KeyError:
                raise e

    @cached_property
    def filenames(self) -> Tuple[Path]:
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

    def _get_max_run_num(self):
        return max(
            int(re.match(S1P.pattern, fl.name).group("run_num"))
            for fl in self.active_contents
        )

    @property
    def max_run_num(self):
        return self._get_max_run_num()

    @classmethod
    def _check_file_consistency(cls, path: Path) -> bool:
        return True

    @classmethod
    def _get_filename_parameters(cls, dct: dict):
        out = {}
        if "repeat_num" not in dct:
            out["repeat_num"] = 1
        return out


class LoadS11(_S11SubDir):
    STANDARD_NAMES = ["Open", "Short", "Match", "External"]
    pattern = r"(?P<load_name>%s)(?P<repeat_num>\d{2})$" % (
        "|".join(LOAD_ALIASES.values())
    )
    known_patterns = (
        f"(?P<load_name>{'|'.join(LOAD_MAPPINGS.keys())})$",
        f"(?P<load_name>{'|'.join(LOAD_ALIASES.values())})$",
        r"(?P<load_name>%s)(?P<repeat_num>\d{1})$" % ("|".join(LOAD_ALIASES.values())),
    )

    known_substitutions = (
        ("AmbientLoad", "Ambient"),
        ("LongCableShort_", "LongCableShorted_"),
    )

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.load_name = LOAD_ALIASES.inverse.get(
            self._match_dict["load_name"], self._match_dict["load_name"]
        )

    def __eq__(self, other):
        return (
            self.__class__.__name__ == other.__class__.__name__
            and self.load_name == other.load_name
        )

    @classmethod
    def _get_filename_parameters(cls, dct: dict):
        out = super()._get_filename_parameters(dct)
        if dct["load_name"] in LOAD_MAPPINGS:
            dct["load_name"] = LOAD_MAPPINGS[dct["load_name"]]
        return out


class AntSimS11(LoadS11):
    pattern = r"(?P<load_name>%s)(?P<repeat_num>\d{2})$" % (
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


class SwitchingState(_S11SubDir):
    pattern = r"(?P<load_name>SwitchingState)(?P<repeat_num>\d{2})$"
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


class ReceiverReading(_S11SubDir):
    pattern = r"(?P<load_name>ReceiverReading)(?P<repeat_num>\d{2})$"
    STANDARD_NAMES = ["Open", "Short", "Match", "ReceiverReading"]
    known_substitutions = (("ReceiverReadings", "ReceiverReading"),)
    known_patterns = ("(?P<load_name>ReceiverReading)",)


class S11Dir(_DataContainer):
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

    def __init__(
        self,
        path: [str, Path],
        *,
        repeat_num: [None, int, dict] = None,
        run_num: [None, int, dict] = None,
        **kwargs,
    ):
        """Class representing the entire S11 subdirectory of an observation

        Parameters
        ----------
        path : str or Path
            Top-level directory of the S11 measurements.
        repeat_num : int or dict, optional
            If int, the repeat num of any applicable sub-directories to use.
            If dict, each key specifies either SwitchingState or ReceiverReading
            and the repeat_num to use for that.
            By default, will find the last repeat.
        run_num : int or dict, optional
            If int, the run num of any applicable sub-directories to use.
            Any given sub-directory uses the same run_num for all files, but each
            sub-directory can use different run_nums
            If dict, each key specifies any of the sub-dir names
            and the repeat_num to use for that.
            By default, will find the last repeat.
        """
        super().__init__(path, **kwargs)

        rep_nums = {}
        for name in (
            ["switching_state", "receiver_reading"]
            + list(LOAD_ALIASES.keys())
            + list(self.get_simulator_names(self.path))
        ):
            try:
                if type(repeat_num) == int:
                    rep_nums[name] = repeat_num
                elif repeat_num is None:
                    rep_nums[name] = self._get_highest_rep_num(
                        self.path, utils.snake_to_camel(name)
                    )
                else:
                    rep_nums[name] = repeat_num.get(
                        name,
                        self._get_highest_rep_num(
                            self.path, utils.snake_to_camel(name)
                        ),
                    )
            except ValueError:
                # Probably no load of that kind in the directory. That's fine.
                pass

        if type(run_num) == int or run_num is None:
            run_nums = {
                **{"switching_state": run_num, "receiver_reading": run_num},
                **{name: run_num for name in LOAD_ALIASES.values()},
            }
        else:
            run_nums = dict(run_num)

        logger.debug(
            f"Highest rep_num for switching state: {self._get_highest_rep_num(self.path, 'SwitchingState')}"
        )

        self.switching_state = SwitchingState(
            self.path / f"SwitchingState{rep_nums['switching_state']:>02}",
            run_num=run_nums.get("switching_state", None),
        )
        self.receiver_reading = ReceiverReading(
            self.path / f"ReceiverReading{rep_nums['receiver_reading']:>02}",
            run_num=run_nums.get("receiver_reading", None),
        )

        for name, load in LOAD_ALIASES.items():
            setattr(
                self,
                name,
                LoadS11(
                    self.path / f"{load}{rep_nums[name]:>02}",
                    run_num=run_nums.get(load, run_nums.get(name, None)),
                ),
            )

        self.simulators = {}
        for name in self.get_simulator_names(path):
            self.simulators[name] = AntSimS11(
                self.path / f"{name}{rep_nums[name]:>02}",
                run_num=run_nums.get(name, None),
            )

    @property
    def load_names(self):
        return tuple(LOAD_ALIASES.keys())

    @property
    def run_num(self):
        """Dictionary specifying run numbers for each load."""
        return {
            k: getattr(self, k).run_num
            for k in list(LOAD_ALIASES.keys()) + ["switching_state", "receiver_reading"]
        }

    @property
    def repeat_num(self):
        return {
            k: getattr(self, k).repeat_num
            for k in list(LOAD_ALIASES.keys()) + ["switching_state", "receiver_reading"]
        }

    @classmethod
    def _get_highest_rep_num(cls, path, kind):
        fls = utils.get_active_files(path)
        fls = [fl for fl in fls if kind in str(fl)]
        rep_nums = [int(str(fl)[-2:]) for fl in fls]
        return max(rep_nums)

    def get_highest_rep_num(self, kind: str):
        """Get the highest repeat number for this kind."""
        return self._get_highest_rep_num(self.path, kind)

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
    def get_simulator_names(cls, path):
        fls = utils.get_active_files(path)
        return {
            fl.name[:-2]
            for fl in fls
            if any(fl.name.startswith(k) for k in ANTENNA_SIMULATORS)
        }

    def __eq__(self, other):
        return self.__class__.__name__ == other.__class__.__name__


class CalibrationObservation(_DataContainer):
    pattern = re.compile(
        r"^Receiver(?P<rcv_num>\d{2})_(?P<temp>\d{2})C_(?P<year>\d{4})_(?P<month>\d{2})_(?P<day>\d{2})_"
        r"(?P<freq_low>\d{3})_to_(?P<freq_hi>\d{3})MHz$"
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
    write_pattern = "Receiver{rcv_num:0>2}_{temp:>02}C_{year:>04}_{month:0>2}_{day:0>2}_{freq_low:0>3}_to_{freq_hi:0>3}MHz"
    _content_type = {
        "S11": S11Dir,
        "Spectra": Spectra,
        "Resistance": Resistances,
        "spectra": Spectra,
        "resistance": Resistances,
        "s11": S11Dir,
    }

    def __init__(
        self,
        path: [str, Path],
        run_num: [int, dict, None] = None,
        repeat_num: [int, None] = None,
        include_previous: bool = True,
        compile_from_def: bool = True,
        **kwargs,
    ):
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
        ambient_temp : int, {15, 25, 35}
            The ambient temperature of the lab measurements.
        run_num : int or dict, optional
            If an integer, the run number to use for all measurements. If None, by default
            uses the last run for each measurement. If a dict, it should specify the
            run number for each applicable measurement.
        repeat_num : int or dict, optional
            If an integer, the repeat number to use for all measurements. If None, by default
            uses the last run for each measurement. If a dict, it should specify the
            repeat number for each applicable measurement.
        fix : bool, optional
            Whether to attempt fixing filenames and the file structure in the observation
            for simple known error cases. Typically it is better to explicitly call
            :method:`check_self` and :method:`check_contents` if you want to perform
            fixes.
        include_previous : bool, optional
            Whether to by default include the previous observation in the same directory
            to supplement the current one if parts are missing.
        compile_from_def : bool, optional
            Whether to attempt compiling a virtual observation from a ``definition.yaml``
            inside the observation directory. This is the default behaviour, but can
            be turned off to enforce that the current directory should be used directly.
        log_level : int, optional
            The file-structure checks can print out a lot of information, which is useful
            when explicitly performing such checks, but typically not desired when one
            simply wants to instantiate the object. Default is to only print errors. Set
            it lower to also print other information.
        """
        # Read the definition file, and combine other observations into a single
        # temporary directory if they exist (otherwise, just symlink this full directory)
        # Note that we need to keep the actual _tmpdir object around otherwise it gets
        # cleaned up!
        self.definition = self.check_definition(Path(path))

        if self.definition.get("entirely_invalid", False):
            logger.warning(
                f"Observation {path} is marked as invalid -- "
                f"proceed with caution! Reason: '{self.definition['entirely_invalid']}'"
            )

        if compile_from_def:
            self._tmpdir, name = self.compile_obs_from_def(path, include_previous)

            path = Path(self._tmpdir.name) / name

        super().__init__(path, **kwargs)

        self._groups = self._match_dict
        self.receiver_num = int(self._groups["rcv_num"])
        self.ambient_temp = int(self._groups["temp"])
        self.year = int(self._groups["year"])
        self.month = int(self._groups["month"])
        self.day = int(self._groups["day"])
        self.freq_low = int(self._groups["freq_low"])
        self.freq_high = int(self._groups["freq_hi"])

        if type(run_num) == int or run_num is None:
            run_nums = {"Spectra": run_num, "Resistance": run_num, "S11": run_num}
        else:
            run_nums = dict(run_num)

        self.spectra = Spectra(
            self.path / "Spectra", run_num=run_nums.get("Spectra", None), **kwargs
        )
        self.resistance = Resistances(
            self.path / "Resistance", run_num=run_nums.get("Resistance", None), **kwargs
        )
        self.s11 = S11Dir(
            self.path / "S11",
            run_num=run_nums.get("S11", None),
            repeat_num=repeat_num,
            **kwargs,
        )

        self.simulator_names = self.get_simulator_names(self.path)

    @classmethod
    def from_observation_yaml(cls, obs_yaml: [str, Path]):
        """Create a CalibrationObservation from a specific YAML format."""
        obs_yaml = Path(obs_yaml)
        assert obs_yaml.exists(), f"{obs_yaml} does not exist!"

        with open(obs_yaml, "r") as fl:
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

        # To keep the temporary directory from being cleaned up, store it on the class.
        cls._tmpdir = tmpdir

        return cls(
            sympath,
            run_num={"S11": s11_run_nums},
            repeat_num=rep_nums,
            fix=False,
            include_previous=False,
            compile_from_def=False,
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

        with open(definition_file, "r") as fl:
            definition = yaml.load(fl, Loader=yaml.FullLoader)

        allowed_keys = {
            "root_obs_dir": str,
            "entirely_invalid": str,
            "include": list,
            "prefer": list,
            "invalid": list,
            "measurements": {
                "resistance_m": {"01": float, "02": float, "03": float},
                "resistance_f": {"01": float, "02": float, "03": float},
            },
            "defaults": {"resistance": dict, "spectra": dict, "s11": dict},
        }

        def _check_grp(defn, allowed):
            for k, v in defn.items():
                if k not in allowed:
                    logger.warning(
                        f"Key {k} found in definitions.yaml, but is not a known keyword."
                    )
                elif isinstance(allowed[k], dict):
                    # Recurse into sub-dictionaries.
                    _check_grp(v, allowed[k])
                elif not isinstance(v, allowed[k]):
                    logger.error(
                        f"Key {k} has wrong type in definitions.yaml. Should be {allowed[k]}, got {type(v)}."
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
    def _validate_match(cls, match: Dict[str, str], filename: str):
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
        """Read all spectra and resistance files. Usually a bad idea."""
        self.spectra.read_all()
        self.resistance.read_all()

    @classmethod
    def get_base_files(cls, path: Path, with_notes=False) -> List[Path]:
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
    ) -> [tempfile.TemporaryDirectory, str]:
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
        ambient_temp : str, optional
            Either '15C', '25C' or '35C'. Actual measurements will be found in a folder
            of this name under ``path``.
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

        # Now make a full symlink directory with these files.
        symdir = tempfile.TemporaryDirectory()

        for fl, fl_abs in files.items():
            sym_path = Path(symdir.name) / obs_name / fl
            if not sym_path.parent.exists():
                sym_path.parent.mkdir(parents=True)
            sym_path.symlink_to(fl_abs)

        return symdir, obs_name

    @classmethod
    def match_path(
        cls, path: Path, root: Path = Path()
    ) -> Tuple[Union[_DataFile, _DataContainer]]:
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
            logger.handlers[0].setLevel(pre_level)

        return parts

    @property
    def run_num(self):
        """Dictionary specifying run numbers for each component"""
        return {
            "S11": self.s11.run_num,
            "Spectra": self.spectra.run_num,
            "Resistance": self.resistance.run_num,
        }

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
