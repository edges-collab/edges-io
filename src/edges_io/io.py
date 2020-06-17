"""
This module defines the overall file structure and internal contents of the
calibration observations. It does *not* implement any algorithms/methods on that data,
making it easier to separate the algorithms from the data checking/reading.
"""

import glob
import logging
import os
import re
import shutil
import tempfile
import warnings
from abc import ABC, abstractmethod
from datetime import datetime
from io import StringIO
from pathlib import Path
from typing import Iterable, List, Tuple, Union

import numpy as np
import read_acq
import toml
import yaml
from bidict import bidict
from cached_property import cached_property

from . import utils
from .data import DATA_PATH
from .h5 import HDF5RawSpectrum
from .logging import logger

LOAD_ALIASES = bidict(
    {
        "ambient": "Ambient",
        "hot_load": "HotLoad",
        "open": "LongCableOpen",
        "short": "LongCableShorted",
    }
)

with open(os.path.join(DATA_PATH, "antenna_simulators.toml")) as fl:
    ANTENNA_SIMULATORS = toml.load(fl)

# Dictionary of misspelled:true mappings.
ANTSIM_REVERSE = {
    v: k for k, val in ANTENNA_SIMULATORS.items() for v in val.get("misspells", [])
}


class _DataFile(ABC):
    def __init__(self, path, fix=False, log_level=40):
        pre_level = logger.getEffectiveLevel()
        logger.setLevel(log_level)

        self.path, match = self.check_self(path)

        logger.setLevel(pre_level)

        try:
            self._match_dict = match.groupdict()
        except AttributeError:
            try:
                self._match_dict = [m.groupdict() for m in match]
            except TypeError:
                self._match_dict = match
        except Exception:
            raise

    @staticmethod
    @abstractmethod
    def check_self(path, fix=False) -> Tuple[Path, Union[re.Match, dict, Iterable]]:
        pass

    @classmethod
    def typestr(cls, name):
        """Generate a string uniquely defining the 'kind of thing' path is.

        The point of this method is to be able to compare two different file/folder
        names to check whether they describe the same kind of thing. For example,
        two Spectrum files from different observations which are both "Ambient" should
        return the same string, even though their dates etc. might be different. However,
        two Spectrum files of different Loads will be different.

        The reason this has to exist (and as a classmethod) is because merely comparing
        the type of the thing is not enough, since the class itself has no knowledge of
        for example the kind of load. But comparing instances is not great either, since
        instances are expected to have a full complement of files to be "valid", but
        one of the main purposes of comparing files is to construct such a full observation.
        """
        return cls.__name__


class _DataContainer(ABC):
    _content_type = None

    def __init__(self, path, fix=False, log_level=40):
        pre_level = logger.getEffectiveLevel()
        logger.setLevel(log_level)

        self.path, match = self.check_self(path, fix)
        self.path = Path(self.path)

        if not self._check_contents_selves(self.path, fix):
            raise utils.FileStructureError()
        if not self._check_all_files_there(self.path):
            raise utils.IncompleteObservation()
        if not self._check_file_consistency(self.path):
            raise utils.InconsistentObservation()

        logger.setLevel(pre_level)

        try:
            self._match_dict = match.groupdict()
        except AttributeError:
            try:
                self._match_dict = [m.groupdict() for m in match]
            except TypeError:
                self._match_dict = match
        except Exception:
            raise

    @classmethod
    @abstractmethod
    def check_self(
        cls, path: Path, fix: bool = False
    ) -> Tuple[Path, Union[re.Match, dict, Iterable]]:
        """
        Check whether the path itself is formatted correctly.

        Parameters
        ----------
        path : Path
            Path to the file/folder.
        definition : dict, optional
            Dictionary specifying the definition of the overall observation, which
            may affect any particular component.
        fix : bool, optional
            Whether to apply straight-forward fixes to the filename format.
        quiet : bool, optional
            Whether to log anything.

        Returns
        -------
        path : Path
            The path to the file, possibly updated by any fix.
        match : dict or None
            If the path was found to be correct, this contains meta information
            contained in the filename.
        """
        pass

    @classmethod
    def check_contents(cls, path: [str, Path], fix=False) -> bool:
        """Abstract method for checking whether the contents of this container are in
         the correct format for the DB"""
        # Check that everything that *is* there has correct format.
        path = Path(path)
        ok_selves = cls._check_contents_selves(path, fix=fix)
        ok_complete = cls._check_all_files_there(path)
        # Check that the files that are there have consistent properties, and are also
        # consistent with outside parameters (eg. if year appears on them, they should
        # be consistent with outer years).
        ok_consistent = cls._check_file_consistency(path)

        return ok_selves and ok_complete and ok_consistent

    @classmethod
    @abstractmethod
    def _check_all_files_there(cls, path: Path) -> bool:
        return True

    @classmethod
    @abstractmethod
    def _check_file_consistency(cls, path: Path) -> bool:
        return True

    @classmethod
    def _check_contents_selves(cls, path: Path, fix=False) -> bool:
        fls = utils.get_active_files(path)

        # Start off with a clean slate for this function.
        logger.errored = 0
        for fl in fls:
            if isinstance(cls._content_type, dict):
                for key, ct in cls._content_type.items():
                    if fl.name.startswith(key):
                        content_type = ct
                        break
                else:
                    logger.error(f"{fl.name} is an extraneous file/folder")

                    if fix:
                        if fl.name == "Notes.odt":
                            shutil.move(fl, fl.with_suffix(".txt"))
                            fl = fl.with_suffix(".txt")
                            logger.success(f"Successfully renamed to {fl}")
                        else:
                            fixed = utils._ask_to_rm(fl)

                            if fixed:
                                logger.success("Successfully removed.")

                    continue
            else:
                content_type = cls._content_type

            fl, _ = content_type.check_self(fl, fix=fix)

            # Recursively check the contents of the contents.
            try:
                content_type.check_contents(fl, fix=fix)
            except AttributeError:
                # It's a DataFile, not a DataContainer
                pass

        ok = not bool(logger.errored)
        logger.errored = 0
        return ok

    @classmethod
    def typestr(cls, name: str) -> str:
        """Generate a string uniquely defining the 'kind of thing' path is.

        The point of this method is to be able to compare two different file/folder
        names to check whether they describe the same kind of thing. For example,
        two Spectrum files from different observations which are both "Ambient" should
        return the same string, even though their dates etc. might be different. However,
        two Spectrum files of different Loads will be different.

        The reason this has to exist (and as a classmethod) is because merely comparing
        the type of the thing is not enough, since the class itself has no knowledge of
        for example the kind of load. But comparing instances is not great either, since
        instances are expected to have a full complement of files to be "valid", but
        one of the main purposes of comparing files is to construct such a full observation.
        """
        return cls.__name__


class _SpectrumOrResistance(_DataFile):
    load_pattern = "|".join(LOAD_ALIASES.values())
    antsim_pattern = "|".join(ANTENNA_SIMULATORS.keys())
    file_pattern = (
        r"(?P<load_name>%s|%s)" % (load_pattern, antsim_pattern)
        + r"_(?P<run_num>\d{2})_(?P<year>\d{4})_(?P<day>\d{3})_("
        r"?P<hour>\d{2})_(?P<minute>\d{2})_(?P<second>\d{2})_lab.(?P<file_format>\w{2,"
        r"3})$"
    )
    supported_formats = []

    def __init__(self, path, fix=False):
        super().__init__(path, fix)

        # Get out metadata
        self._groups = self._match_dict

    @classmethod
    def typestr(cls, name: str):
        return cls.__name__ + re.match(cls.file_pattern, name).groupdict()["load_name"]

    @classmethod
    def _fix(cls, root, basename):
        if "AmbientLoad" in basename:
            newname = basename.replace("AmbientLoad", "Ambient")
        elif "LongCableShort_" in basename:
            newname = basename.replace("LongCableShort_", "LongCableShorted_")
        else:
            newname = basename

        match = re.search(cls.file_pattern, newname)

        if match is not None:
            shutil.move(os.path.join(root, basename), os.path.join(root, newname))
            logger.success("Successfully converted to {}".format(newname))
            return os.path.join(root, newname), match

        bad_patterns = [
            (
                r"^(?P<load_name>%s|%s)" % (cls.load_pattern, cls.antsim_pattern)
                + r"_25C_(?P<month>\d{1,2})_(?P<day>\d{1,2})_("
                r"?P<year>\d\d\d\d)_(?P<hour>\d{1,2})_(?P<minute>\d{"
                r"1,2})_(?P<second>\d{1,2}).(?P<file_format>\w{2,3})$"
            ),
            (
                "^(?P<load_name>{})".format(cls.load_pattern)
                + r"(?P<run_num>\d{1,2})_25C_(?P<month>\d{1,"
                r"2})_(?P<day>\d{1,2})_(?P<year>\d\d\d\d)_("
                r"?P<hour>\d{1,2})_(?P<minute>\d{1,"
                r"2})_(?P<second>\d{1,2}).(?P<file_format>\w{2,3})$"
            ),
            (
                r"(?P<load_name>%s|%s)" % (cls.load_pattern, cls.antsim_pattern)
                + r"_(?P<year>\d{4})_(?P<day>\d{3})_"
                r"(?P<hour>\d{2})_(?P<minute>\d{2})_(?P<second>\d{2})_lab."
                r"(?P<file_format>\w{2,3})$"
            ),
            (
                r"(?P<load_name>%s|%s)" % (cls.load_pattern, cls.antsim_pattern)
                + r"_(?P<year>\d{4})_(?P<day>\d{3})_(?P<hour>\d{2}).(?P<file_format>\w{2,3})$"
            ),
            (
                r"(?P<load_name>%s|%s)" % (cls.load_pattern, cls.antsim_pattern)
                + r"_(?P<year>\d{4})_(?P<day>\d{3})_(?P<hour>\d{2}).(?P<file_format>\w{2,3})$"
            ),
            (
                r"(?P<load_name>%s|%s)" % (cls.load_pattern, cls.antsim_pattern)
                + r"_(?P<year>\d{4})_(?P<day>\d{3})_lab.(?P<file_format>\w{2,3})$"
            ),
            (
                r"(?P<load_name>%s|%s)" % (cls.load_pattern, cls.antsim_pattern)
                + r"_(?P<run_num>\d)_(?P<year>\d{4})_(?P<day>\d{3})_lab.(?P<file_format>\w{2,"
                r"3})$"
            ),
            (
                r"(?P<load_name>%s|%s)" % (cls.load_pattern, cls.antsim_pattern)
                + r"_\d{2}C_(?P<month>\d{1,2})_(?P<day>\d{1,2})_(?P<year>\d{4})_(?P<hour>\d{"
                r"1,2})_(?P<minute>\d{1,2})_(?P<second>\d{1,2}).(?P<file_format>\w{2,3})"
            ),
            (
                r"(?P<load_name>%s)" % ("|".join(ANTSIM_REVERSE.keys()))
                + r"_\d{2}C_(?P<month>\d{1,2})_(?P<day>\d{1,2})_(?P<year>\d{4})_(?P<hour>\d{"
                r"1,2})_(?P<minute>\d{1,2})_(?P<second>\d{1,2}).(?P<file_format>\w{2,3})"
            ),
        ]
        i = 0
        while match is None and i < len(bad_patterns):
            pattern = bad_patterns[i]
            match = re.search(pattern, newname)
            i += 1

        if match is None:
            logger.warning("\tCould not auto-fix it.")

            fixed = utils._ask_to_rm(os.path.join(root, basename))
            if fixed:
                logger.success("Successfully removed.")
            return None, None
        else:
            dct = match.groupdict()

            if "month" in dct:
                jd = utils.ymd_to_jd(
                    match.group("year"), match.group("month"), match.group("day")
                )
            else:
                jd = dct["day"]

            if "run_num" not in dct:
                dct["run_num"] = "01"

            if "hour" not in dct:
                dct["hour"] = "00"

            if "minute" not in dct:
                dct["minute"] = "00"

            if "second" not in dct:
                dct["second"] = "00"

            # Switch Antenna Simulator "misspells" to true form.
            if dct["load_name"] in ANTSIM_REVERSE:
                dct["load_name"] = ANTSIM_REVERSE[dct["load_name"]]

            newname = (
                "{load_name}_{run_num:0>2}_{year:0>4}_{jd:0>3}_{hour:0>2}_{minute:0>2}_"
                "{second:0>2}_lab.{file_format}".format(jd=jd, **dct)
            )
            newpath = os.path.join(root, newname)

            match = re.search(cls.file_pattern, newname)

            if match is not None:
                logger.success(f"Successfully converted to {newname}")
                shutil.move(os.path.join(root, basename), newpath)
                return newpath, match
            else:
                return None, None

    @classmethod
    def check_self(cls, path, fix=False):
        if isinstance(path, (str, Path)):
            path = [Path(path)]
        else:
            path = [Path(p) for p in path]

        base_fnames = [fname.name for fname in path]
        root = path[0].parent

        matches = []
        for i, basename in enumerate(base_fnames):
            match = re.search(cls.file_pattern, basename)
            if match is None:
                logger.error(
                    f"The file {basename} does not have the correct format for a {cls.__name__}"
                )

                if fix:
                    newname, match = cls._fix(root, basename)
                    path[i] = newname

            if match is not None:
                newname = os.path.basename(path[i])
                groups = match.groupdict()
                if int(groups["run_num"]) < 1:
                    logger.error(f"The run_num for {newname} is less than one!")
                if not (2010 <= int(groups["year"]) <= 2030):
                    logger.error(
                        f"The year for {newname} ({groups['year']}) is a bit strange!"
                    )
                if not (0 <= int(groups["day"]) <= 366):
                    logger.error(
                        f"The day for {newname} ({groups['day']}) is outside the number "
                        f"of days in a year"
                    )
                if not (0 <= int(groups["hour"]) <= 24):
                    logger.error(f"The hour for {newname} is outside 0-24!")
                if not (0 <= int(groups["minute"]) <= 60):
                    logger.error(f"The minute for {newname} is outside 0-60!")
                if not (0 <= int(groups["second"]) <= 60):
                    logger.error(f"The second for {newname} is outside 0-60!")
                if groups["file_format"] not in cls.supported_formats:
                    logger.error(
                        f"The file {newname} is not of a supported format "
                        f"({cls.supported_formats}). Got format {groups['file_format']}"
                    )
                matches.append(match)
        return path, matches

    @classmethod
    def from_load(cls, load, direc, run_num=None, filetype=None):
        """
        Initialize the object in a simple way.

        Parameters
        ----------
        load : str
            The load name (eg. 'Ambient', 'HotLoad') or its alias (eg. 'ambient',
            'hot_load').
        direc : str
            The directory in which to search for relevant data
        run_num : int, optional
            The run number of the data to use. Default, the last run. Each run is
            independent and different run_nums may be used for different loads.
        filetype : str, optional
            The filetype of the data. Must be one of the supported formats. Defaults
            to `_default_filetype`.
        """
        if load in LOAD_ALIASES:
            load = LOAD_ALIASES[load]

        if load not in LOAD_ALIASES.values() and load not in ANTENNA_SIMULATORS:
            logger.error(
                "The load specified {} is not one of the options available.".format(
                    load
                )
            )

        files = glob.glob(
            os.path.join(direc, "{load}_??_????_???_??_??_??_lab.*".format(load=load))
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
                restricted_files = [fl for fl in files if fl.endswith("." + ftype)]
                if restricted_files:
                    break
            files = restricted_files

            if not files:
                raise ValueError(
                    f"No files exist for the load {load} for any filetype on that path: {direc}"
                )

        # Restrict to the given run_num (default last run)
        run_nums = [
            int(os.path.basename(fl)[len(load) + 1 : len(load) + 3]) for fl in files
        ]
        if run_num is None:
            run_num = max(run_nums)

        pre_files = files.copy()
        files = [fl for fl, num in zip(files, run_nums) if num == run_num]

        if not files:
            raise ValueError(
                "No {} files exist on path ({}) with run_num={}. Potential files: {}".format(
                    load, direc, run_num, pre_files
                )
            )

        return cls(files)

    @cached_property
    def run_num(self):
        """The run number of the data. All run_nums must be the same for all files in
        the data.

        Every observation may have several runs. Note that different runs may be mixed
        for different loads.
        """
        # Ensure all load names are the same
        if any(
            group["run_num"] != self._groups[0]["run_num"] for group in self._groups
        ):
            raise IOError("Two files given with incompatible run_nums")
        return self._groups[0]["run_num"]

    @cached_property
    def year(self):
        """Year on which data acquisition began"""
        # Ensure all load names are the same
        if any(group["year"] != self._groups[0]["year"] for group in self._groups):
            raise IOError("Two files given with incompatible years")
        return int(self._groups[0]["year"])

    @cached_property
    def days(self):
        """List of integer days (one per file) at which data acquisition was begun"""
        days = [int(group["day"]) for group in self._groups]
        if max(days) - min(days) > 30:
            logger.warning(
                "Spectra taken suspiciously far apart [{} days]".format(
                    max(days) - min(days)
                )
            )

        return days

    @cached_property
    def load_name(self):
        return LOAD_ALIASES.inverse.get(
            self._groups[0]["load_name"], self._groups[0]["load_name"]
        )

    @cached_property
    def hours(self):
        """List of integer hours (one per file) at which data acquisition was begun"""
        return [int(group["hour"]) for group in self._groups]

    @cached_property
    def minutes(self):
        """List of integer minutes (one per file) at which data acquisition was begun"""
        return [int(group["minute"]) for group in self._groups]

    @cached_property
    def seconds(self):
        """List of integer seconds (one per file) at which data acquisition was begun"""
        return [int(group["second"]) for group in self._groups]

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
    def file_format(self):
        """The file format of the data to be read."""
        formats = [os.path.splitext(fl)[1][1:] for fl in self.path]
        if any(
            format != formats[0] or format not in self.supported_formats
            for format in formats
        ):
            raise ValueError("not all file formats are the same!")
        return formats[0]

    @cached_property
    def data(self) -> [HDF5RawSpectrum, List[HDF5RawSpectrum]]:
        """A view of the data in the file as a HDF5Object.

        If the file is an ACQ file, it will be read completely into memory and cast
        into the same format as a :class:`~h5.HDF5Object` so that the API is the same.

        If the number of files is more than one, `data` will be a list of objects.
        """
        objs = []
        for fl in self.path:
            if self.file_format == "h5":
                objs.append(HDF5RawSpectrum(fl))
            elif self.file_format == "acq":
                spec, freq, time, meta = self._read_acq(fl)
                objs.append(
                    HDF5RawSpectrum.from_data(
                        {
                            "spectra": spec,
                            "time_ancillary": time,
                            "freq_ancillary": freq,
                            "meta": meta,
                        }
                    )
                )

        if len(objs) == 1:
            return objs[0]
        else:
            return objs

    #
    # def read(self):
    #     """
    #     Read the files of the object, and concatenate their data.
    #
    #     .. warning:: Deprecated.
    #
    #     Adds the attributes 'p0', 'p1', 'p2' and 'Qratio'.
    #     """
    #     warnings.warn(
    #         "Do not use this function, it will soon be removed. "
    #         "Instead access data through the 'data' attribute",
    #         category=DeprecationWarning
    #     )
    #
    #     out = {}
    #     keys = ["p0", "p1", "p2", "Qratio"]
    #     for fl in self.path:
    #         spectra, freq, tm, meta = getattr(self, "_read_" + self.file_format)(fl)
    #
    #         for key in keys:
    #             if key not in out:
    #                 out[key] = spectra[key]
    #             else:
    #                 out[key] = np.concatenate((out[key], spectra[key]), axis=1)
    #     setattr(self, "spectra", out)
    #     setattr(self, "metadata", meta)
    #     setattr(self, 'time_ancillary', tm)
    #     setattr(self, 'freq_ancillary', freq)
    #     return out, freq, tm, meta

    @staticmethod
    def _read_acq(file_name):
        Q, px, anc = read_acq.decode_file(
            file_name, progress=False, write_formats=[], meta=True
        )

        freq_anc = {"frequencies": anc.frequencies}
        time_anc = {
            name: anc["time_data"][name] for name in anc["time_data"].dtype.names
        }
        spectra = {"Qratio": Q.T, "p0": px[0].T, "p1": px[1].T, "p2": px[2].T}

        meta = anc.meta
        return spectra, freq_anc, time_anc, meta

    # @staticmethod
    # def _read_h5(file_name):
    #     spectra = {}
    #     freq_anc = {}
    #     time_anc = {}
    #
    #     with h5py.File(file_name, "r") as fl:
    #         spectra["Qratio"] = fl["Qratio"][...]
    #         spectra["p0"] = fl["p0"][...]
    #         spectra["p1"] = fl["p1"][...]
    #         spectra["p2"] = fl["p2"][...]
    #
    #         meta = dict(fl.attrs)
    #
    #         for k, v in fl['freq_ancillary'].items():
    #             freq_anc[k] = v[...]
    #
    #         for k, v in fl['time_ancillary'].items():
    #             time_anc[k] = v[...]
    #
    #     return spectra, freq_anc, time_anc, meta


class Resistance(_SpectrumOrResistance):
    """An object representing a resistance measurement (and its structure)."""

    supported_formats = ("csv",)

    def __init__(self, path, fix=False, store_data=True):
        super(Resistance, self).__init__(path, fix=fix)
        self.path = self.path[0]
        self.store_data = store_data

    @classmethod
    def check_self(cls, path, fix=False):
        path, match = super(Resistance, cls).check_self(path, fix)
        if len(path) > 1:
            logger.error(
                "Only one resistance file should exist for each load and run_num"
            )
        return path, match

    @cached_property
    def file_format(self):
        """The file format of the data to be read."""
        return "csv"

    def read(self):

        try:
            return self._data, self._meta
        except AttributeError:
            with open(self.path, "r", errors="ignore") as fl:
                if fl.readline().startswith("FLUKE"):
                    data, meta = self.read_old_style_csv(self.path)
                else:
                    data, meta = self.read_new_style_csv(self.path)

            if self.store_data:
                self._data = data
                self._meta = meta

            return self._data, self._meta

    @classmethod
    def read_new_style_csv(cls, path: [str, Path]):
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
    def read_old_style_csv(cls, path):
        # Weirdly, some old-style files use KOhm, and some just use Ohm.

        # These files have bad encoding, which we can ignore. This means we have to
        # read in the whole thing as text first (while ignoring errors) and construct
        # a StringIO object to pass to genfromtxt.
        with open(path, "r", errors="ignore") as fl:
            while not fl.readline().startswith("Start Time,"):
                continue

            # Get the total number of actual lines
            nlines = int(fl.readline().split(",")[4])

            # The last line can be only half there, so omit it
            nlines -= 1

            while not fl.readline().startswith("Reading,"):
                continue

            s = StringIO("".join([next(fl) for i in range(nlines)]))

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


class _SpectraOrResistanceFolder(_DataContainer):
    folder_pattern = None

    def __init__(self, path, run_num=None, filetype=None, fix=False):
        """Collection of spectra in an observation"""
        super().__init__(path, fix)

        if type(run_num) is int or run_num is None:
            run_nums = {load: run_num for load in LOAD_ALIASES.values()}
        else:
            run_nums = run_num

        for name, load in LOAD_ALIASES.items():
            setattr(
                self,
                name,
                self._content_type.from_load(
                    load, path, run_nums.get(load, None), filetype
                ),
            )

        # Populate simulators.
        self.simulators = {}
        for name in self.get_simulator_names(self.path):
            self.simulators[name] = self._content_type.from_load(
                name, self.path, run_nums.get(name, None), filetype
            )

    @classmethod
    def check_self(cls, path, fix=False):
        logger.structure("Checking {} folder contents at {}".format(cls.__name__, path))

        match = re.search(cls.folder_pattern, os.path.basename(path))
        if match is None:
            logger.error(
                "{} directory should be called {}".format(
                    cls.__name__, cls.folder_pattern
                )
            )

        return path, match

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
            re.search(cls._content_type.file_pattern, fl.name).groupdict() for fl in fls
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
    folder_pattern = "Spectra"
    _content_type = Spectrum


class Resistances(_SpectraOrResistanceFolder):
    folder_pattern = "Resistance"
    _content_type = Resistance


class S1P(_DataFile):
    file_pattern = r"(?P<kind>\w+)(?P<run_num>\d{2}).s1p$"

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

    @classmethod
    def typestr(cls, name: str) -> str:
        return cls.__name__ + re.match(cls.file_pattern, name).groupdict()["kind"]

    def __init__(self, path, fix=False):
        super().__init__(path, fix)

        self.kind = self._match_dict["kind"]
        self.run_num = int(self._match_dict["run_num"])

        self.s11, self.freq = self.read(self.path)

    @classmethod
    def check_self(cls, path, fix=False):
        basename = os.path.basename(path)
        match = re.search(cls.file_pattern, basename)
        if match is None:
            logger.error(
                "The file {} has the wrong filename format for an S11 file".format(path)
            )
        else:
            groups = match.groupdict()
            if groups["kind"] not in cls.POSSIBLE_KINDS:
                logger.error(
                    "The file {} has a kind ({}) that is not supported. "
                    "Possible: {}".format(path, groups["kind"], cls.POSSIBLE_KINDS)
                )

                if fix and groups["kind"].capitalize() in cls.POSSIBLE_KINDS:
                    shutil.move(
                        path, os.path.join(os.path.dirname(path), basename.capitalize())
                    )
                    logger.success("Changed to ", basename.capitalize())

            if int(groups["run_num"]) < 1:
                logger.error(
                    "The file {} has a run_num ({}) less than one".format(
                        path, groups["run_num"]
                    )
                )
        return path, match

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
    folder_pattern = None

    @classmethod
    def typestr(cls, name: str) -> str:
        return (
            cls.__name__ + re.match(cls.folder_pattern, name).groupdict()["load_name"]
        )

    def __init__(self, path, run_num=None, fix=False):
        super().__init__(path, fix)

        self.run_num = run_num or self._get_max_run_num()

        for name in self.STANDARD_NAMES:
            setattr(
                self,
                name.lower(),
                S1P(os.path.join(path, name + "{:>02}.s1p".format(self.run_num))),
            )

        # All frequencies should be the same.
        self.freq = getattr(self, self.STANDARD_NAMES[0].lower()).freq

        self.filenames = [
            getattr(self, thing.lower()).path for thing in self.STANDARD_NAMES
        ]

    @property
    def active_contents(self):
        return utils.get_active_files(self.path)

    @classmethod
    def check_self(cls, path, fix=False):
        if not os.path.exists(path):
            logger.error("The path {} does not exist!".format(path))
            return path, None

        match = re.search(cls.folder_pattern, os.path.basename(path))

        if match is None:
            logger.error(
                "The folder {} did not match any of the correct folder name "
                "criteria. Required pattern: {} [class={}]".format(
                    path, cls.folder_pattern, cls.__name__
                )
            )

            if fix:
                if "AmbientLoad" in path:
                    newpath = path.replace("AmbientLoad", "Ambient")
                elif "LongCableShort_" in path or path.endswith("LongCableShort"):
                    newpath = path.replace("LongCableShort", "LongCableShorted")
                elif "InternalSwitch" in path:
                    newpath = path.replace("InternalSwitch", "SwitchingState")
                else:
                    newpath = path

                # Sometimes they don't have repeat_nums attached.
                if os.path.basename(newpath) == "SwitchingState":
                    newpath = newpath.replace("SwitchingState", "SwitchingState01")
                elif os.path.basename(newpath) == "ReceiverReading":
                    newpath = newpath.replace("ReceiverReading", "ReceiverReading01")

                if newpath != path:
                    shutil.move(path, newpath)
                    path = newpath

                match = re.search(cls.folder_pattern, os.path.basename(path))

                if match is not None:
                    logger.success("Successfully converted to {}".format(path))

        return path, match

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
            int(re.match(S1P.file_pattern, fl.name).group("run_num"))
            for fl in self.active_contents
        )

    @classmethod
    def _check_file_consistency(cls, path: Path) -> bool:
        return True


class LoadS11(_S11SubDir):
    STANDARD_NAMES = ["Open", "Short", "Match", "External"]
    folder_pattern = "(?P<load_name>{})$".format("|".join(LOAD_ALIASES.values()))

    def __init__(self, direc, run_num=None, fix=False):
        super().__init__(direc, run_num, fix)
        self.load_name = LOAD_ALIASES.inverse.get(
            self._match_dict["load_name"], self._match_dict["load_name"]
        )

    def __eq__(self, other):
        return (
            self.__class__.__name__ == other.__class__.__name__
            and self.load_name == other.load_name
        )


class AntSimS11(LoadS11):
    folder_pattern = r"(?P<load_name>%s)$" % ("|".join(ANTENNA_SIMULATORS.keys()))

    @classmethod
    def check_self(cls, path, fix=False):
        path, match = super().check_self(path, fix)

        if match is None and fix:
            bad_patterns = [r"(?P<load_name>%s)$" % ("|".join(ANTSIM_REVERSE.keys()))]

            for pattern in bad_patterns:
                match = re.search(pattern, os.path.basename(path))

            if match is not None:
                loadname = match.groupdict()["load_name"]

                if loadname in ANTSIM_REVERSE:
                    loadname = ANTSIM_REVERSE[loadname]

                newpath = os.path.join(os.path.dirname(path), loadname)

                shutil.move(path, newpath)
                path = newpath
                logger.success("Successfully converted to {}".format(path))

        return path, match


class _RepeatNumberableS11SubDir(_S11SubDir):
    def __init__(self, direc, run_num=None, fix=False):
        super().__init__(direc, run_num, fix)
        self.repeat_num = int(self._match_dict["repeat_num"])

    def __eq__(self, other):
        return self.__class__.__name__ == other.__class__.__name__

    @classmethod
    def typestr(cls, name: str) -> str:
        return cls.__name__


class SwitchingState(_RepeatNumberableS11SubDir):
    folder_pattern = r"SwitchingState(?P<repeat_num>\d{2})$"
    STANDARD_NAMES = [
        "Open",
        "Short",
        "Match",
        "ExternalOpen",
        "ExternalShort",
        "ExternalMatch",
    ]


class ReceiverReading(_RepeatNumberableS11SubDir):
    folder_pattern = r"ReceiverReading(?P<repeat_num>\d{2})$"
    STANDARD_NAMES = ["Open", "Short", "Match", "ReceiverReading"]


class S11Dir(_DataContainer):
    _content_type = {
        **{load: LoadS11 for load in LOAD_ALIASES.values()},
        **{
            "SwitchingState": SwitchingState,
            "ReceiverReading": ReceiverReading,
            "InternalSwitch": SwitchingState,  # To catch the old way so it can be fixed.
            "LongCableShort": LoadS11,
        },
        **{key: AntSimS11 for key in ANTENNA_SIMULATORS.keys()},
        **{key: AntSimS11 for key in ANTSIM_REVERSE.keys()},
    }

    def __init__(self, path: [str, Path], repeat_num=None, run_num=None, fix=False):
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
        super().__init__(path, fix)

        if type(repeat_num) == int:
            sw_rep_num = repeat_num
            rr_rep_num = repeat_num
        elif repeat_num is None:
            sw_rep_num = self._get_highest_rep_num(path, "SwitchingState")
            rr_rep_num = self._get_highest_rep_num(path, "ReceiverReading")
        else:
            sw_rep_num = repeat_num["SwitchingState"]
            rr_rep_num = repeat_num["ReceiverReading"]

        if type(run_num) == int or run_num is None:
            run_nums = {
                **{"SwitchingState": run_num, "ReceiverReading": run_num},
                **{name: run_num for name in LOAD_ALIASES.values()},
            }
        else:
            run_nums = run_num

        logger.debug(
            "Highest rep_num for switching state: {}".format(
                self._get_highest_rep_num(path, "SwitchingState")
            )
        )

        self.switching_state = SwitchingState(
            os.path.join(path, "SwitchingState{:>02}".format(sw_rep_num)),
            run_num=run_nums.get("SwitchingState", None),
        )
        self.receiver_reading = ReceiverReading(
            os.path.join(path, "ReceiverReading{:>02}".format(rr_rep_num)),
            run_num=run_nums.get("ReceiverReading", None),
        )

        for name, load in LOAD_ALIASES.items():
            setattr(
                self,
                name,
                LoadS11(os.path.join(path, load), run_num=run_nums.get(load, None)),
            )

        self.simulators = {}
        for name in self.get_simulator_names(path):
            self.simulators[name] = AntSimS11(
                os.path.join(path, name), run_num=run_nums.get(name, None)
            )

    @classmethod
    def _get_highest_rep_num(cls, path, kind):
        fls = utils.get_active_files(path)
        fls = [fl for fl in fls if kind in str(fl)]
        rep_nums = [int(str(fl)[-2:]) for fl in fls]
        return max(rep_nums)

    @classmethod
    def check_self(cls, path, fix=False):
        logger.structure("Checking S11 folder contents at {}".format(path))

        if not os.path.exists(path):
            logger.error(f"This path does not exist: {path}")

        if os.path.basename(path) != "S11":
            logger.error("The S11 folder should be called S11")

        return path, os.path.basename(path) == "S11"

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
            fl.name
            for fl in fls
            if any(fl.name.startswith(k) for k in ANTENNA_SIMULATORS)
        }

    def __eq__(self, other):
        return self.__class__.__name__ == other.__class__.__name__


class CalibrationObservation(_DataContainer):
    file_pattern = re.compile(
        r"^Receiver(?P<rcv_num>\d{2})_(?P<year>\d{4})_(?P<month>\d{2})_(?P<day>\d{2})_("
        r"?P<freq_low>\d{3})_to_(?P<freq_hi>\d{3})MHz/(?P<temp>\d{2})C$"
    )
    _content_type = {"S11": S11Dir, "Spectra": Spectra, "Resistance": Resistances}

    def __init__(
        self,
        path: [str, Path],
        ambient_temp: int = 25,
        run_num: [int, dict, None] = None,
        repeat_num: [int, None] = None,
        fix: bool = False,
        include_previous: bool = True,
        compile_from_def: bool = True,
        log_level=40,
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
        if ambient_temp not in [15, 25, 35]:
            raise ValueError("ambient temp must be one of 15, 25, 35!")

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
            self._tmpdir, name = self.compile_obs_from_def(
                path, f"{ambient_temp}C", include_previous
            )

            path = Path(self._tmpdir.name) / name / f"{ambient_temp}C"
        else:
            path = Path(path) / f"{ambient_temp}C"

        super().__init__(path, fix)

        self.ambient_temp = ambient_temp
        self._groups = self._match_dict
        self.receiver_num = int(self._groups["rcv_num"])
        self.year = int(self._groups["year"])
        self.month = int(self._groups["month"])
        self.day = int(self._groups["day"])
        self.freq_low = int(self._groups["freq_low"])
        self.freq_high = int(self._groups["freq_hi"])

        if type(run_num) == int or run_num is None:
            run_nums = {"Spectra": run_num, "Resistance": run_num, "S11": run_num}
        else:
            run_nums = run_num

        self.spectra = Spectra(
            self.path / "Spectra", run_num=run_nums.get("Spectra", None), fix=fix
        )
        self.resistance = Resistances(
            self.path / "Resistance", run_num=run_nums.get("Resistance", None), fix=fix
        )
        self.s11 = S11Dir(
            self.path / "S11",
            run_num=run_nums.get("S11", None),
            repeat_num=repeat_num,
            fix=fix,
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

        sympath = (
            Path(tmpdir.name)
            / f"Receiver{meta['receiver']:02d}_{meta['year']}_{meta['month']:02d}_{meta['day']:02d}_040_to_200MHz/25C"
        )
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

        kind_map = {
            "ambient": "Ambient",
            "hot_load": "HotLoad",
            "open": "LongCableOpen",
            "short": "LongCableShorted",
            "receiver": "ReceiverReading01",
            "switch": "SwitchingState01",
        }

        # Symlink the S11 files.
        for key, val in files["s11"].items():
            direc = s11 / kind_map.get(key, utils.snake_to_camel(key))
            direc.mkdir()

            for key2, val2 in val.items():
                if key2 == "receiver":
                    key2 = "receiver_reading"

                filename = utils.snake_to_camel(key2) + "01.s1p"
                (direc / filename).symlink_to(root / val2)

        # To keep the temporary directory from being cleaned up, store it on the class.
        cls._tmpdir = tmpdir

        return cls(
            sympath.parent,
            ambient_temp=25,
            run_num=1,
            repeat_num=1,
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

                if key == "s11":
                    for key3 in ["open", "short", "match", "external"]:
                        assert (
                            key3 in files[key][key2]
                        ), f"{key3} must be in observation YAML 'files.{key}.{key2}'"
                else:
                    for fl in files[key][key2]:
                        assert (
                            len(list(root.glob(fl))) > 0
                        ), f"File '{fl}' included at files.{key}.{key2} does not exist or match any glob patterns."

            if key == "s11":
                for key2 in ["receiver", "switch"]:
                    assert (
                        key2 in files[key]
                    ), f"{key2} must be in observation YAML 'files.{key}'. Available: {list(files[key].keys())}"

                    for key3 in ["match", "open", "short"]:
                        assert (
                            key3 in files[key][key2]
                        ), f"{key3} must be in observation YAML 'files.{key}.{key2}'"
                        assert (
                            root / files[key][key2][key3]
                        ).exists(), f"File '{files[key][key2][key3]}' included at files.{key}.{key2}.{key3} does not exist."

                    if key2 == "receiver":
                        assert (
                            "receiver" in files[key][key2]
                        ), f"'receiver' must be in observation YAML 'files.{key}.{key2}'"
                        assert (
                            root / files[key][key2]["receiver"]
                        ).exists(), f"File {files[key][key2]['receiver']} included at files.{key}.{key2}.receiver does not exist."
                    elif key2 == "switch":
                        for key3 in [
                            "external_match",
                            "external_open",
                            "external_short",
                        ]:
                            assert (
                                key3 in files[key][key2]
                            ), f"{key3} must be in observation YAML 'files.{key}.{key2}'"
                            assert (root / files[key][key2][key3]).exists(), (
                                f"File '{files[key][key2][key3]}' included at "
                                f"files.{key}.{key2}.{key3} does not exist."
                            )

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
            "measurements": {"resistance_m": float, "resistance_f": float},
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
    def check_self(cls, path: [str, Path], fix: bool = False):
        path = Path(path).absolute()

        logger.structure(f"Checking root folder: {path}")

        if not path.exists():
            raise IOError(f"The path {path} does not exist!")

        base = path.parents[1]  # Gets root obs dir.
        name = str(Path(path.parent.name) / Path(path.name))

        # Warn if this is an invalid observation entirely. Also, we don't check the
        # observation then, as it's annoyingly difficult.
        if path.parent.suffix in [".invalid", ".old"]:
            logger.warning(
                f"Observation {path.parent.name} is marked as invalid -- "
                f"proceed with caution!"
            )
            return path, None

        # If it was specified as entirely valid in the yaml, we can continue to check
        # other aspects, but we'll still emit a warning

        match = cls.file_pattern.search(name)

        if match is None:
            logger.error(
                f"Calibration Observation directory name is in the wrong format! Got {name}"
            )

            if fix:

                bad_pattern = re.compile(
                    r"^Receiver(\d{1,2})_(\d{4})_(\d{1,2})_(\d{1,2})_(\d{2,3})_to_(\d{"
                    r"2,3})MHz/(?P<temp>\d{2})C$"
                )

                match = bad_pattern.search(name)

                if match is None:
                    bad_pattern = re.compile(
                        r"Receiver(?P<rcv_num>\d{2})_(?P<year>\d{4})_(?P<month>\d{2})_"
                        r"(?P<day>\d{2})_(?P<freq_low>\d{3})"
                        r"_to_(?P<freq_hi>\d{3})_MHz/(?P<temp>\d{2})C$"
                    )
                    match = bad_pattern.search(name)

                if match is not None:
                    new_name = "Receiver{:0>2}_{}_{:0>2}_{:0>2}_{:0>3}_to_{:0>3}MHz/{}C".format(
                        *match.groups()
                    )
                    shutil.move(
                        os.path.normpath(os.path.dirname(path)),
                        os.path.join(base, os.path.dirname(new_name)),
                    )

                    # # If top-level directory is now empty, remove it.
                    # if not glob.glob(os.path.join(os.path.dirname(os.path.normpath(path))), "*"):
                    #     os.rmdir(os.path.dirname(os.path.normpath(path)))

                    name = new_name
                    path = base / name
                    logger.success(f"Successfully renamed to {new_name}")
                else:
                    logger.warning("Failed to fix the name scheme")

        if match is not None:
            groups = match.groupdict()
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

            logger.info("Calibration Observation Metadata:")
            for k, v in groups.items():
                logger.info(f"\t{k}: {v}")

        logger.debug(f"\tReturning path={path}")

        return path, match

    @classmethod
    def path_to_datetime(cls, path: [str, Path]):
        pre_level = logger.getEffectiveLevel()
        logger.setLevel(39)
        try:
            path, match = cls.check_self(path, fix=False)
        except Exception as e:
            raise (e)
        finally:
            logger.setLevel(pre_level)

        if match:
            grp = match.groupdict()
            return datetime(int(grp["year"]), int(grp["month"]), int(grp["day"]))
        else:
            raise utils.FileStructureError("The path is not valid for an Observation.")

    @classmethod
    def _check_all_files_there(cls, path: Path) -> bool:
        ok = True
        for folder in ["S11", "Spectra", "Resistance"]:
            if not (path / folder).exists():
                logger.error(f"No {folder} folder in observation!")
                ok = False
        return ok

    @classmethod
    def _check_file_consistency(cls, path: Path) -> bool:
        # checks whether simulators are the same between each.
        cls.get_simulator_names(path)
        return True

    @classmethod
    def get_simulator_names(cls, path):
        # Go through the subdirectories and check their simulators
        dct = {
            name: tuple(sorted(kls.get_simulator_names(os.path.join(path, name))))
            for name, kls in cls._content_type.items()
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
        files = utils.get_file_list(
            path,
            filter=lambda x: x.suffix not in [".invalid", ".old"]
            and x.name not in other_ignores,
            ignore=invalid,
        )

        return files

    @classmethod
    def compile_obs_from_def(
        cls, path: Path, ambient_temp="25C", include_previous=True
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
        path = Path(path)
        obs_name = path.name
        path = (path / ambient_temp).absolute()

        assert path.exists(), f"{path} does not exist!"
        definition = cls.check_definition(path)

        # Now include files from other observations if they don't already exist.
        root_obs = definition.get("root_obs_dir", None)

        if root_obs is None:
            root_obs = path.parents[1]
        else:
            # Root observation directory should be relative to the definition file.
            if not Path(root_obs).is_absolute():
                root_obs = (path / root_obs).resolve()

        files = {fl.relative_to(path.parents[1]): fl for fl in cls.get_base_files(path)}

        file_parts = {
            fl.relative_to(obs_name): cls.match_path(fl, root=path.parents[1])
            for fl in files
        }
        # Actually need files to *not* have the top-level name
        files = {fl.relative_to(obs_name): fl_abs for fl, fl_abs in files.items()}

        def _include_extra(roots, prefer):
            for inc_path in roots:
                # Need to get this root_obs if inc_path is absolute, because we need
                # to know where the observation starts (in the path)
                if inc_path.is_absolute():
                    indx = inc_path.parts.index(ambient_temp)
                    this_root_obs = inc_path.parents[len(inc_path.parts) - indx]
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
        path_parts[0] = os.path.join(path_parts[0], path_parts[1])
        del path_parts[1]

        try:
            parts = tuple()
            full_part = root
            for part in path_parts:
                full_part = Path(full_part) / Path(part)

                for thing in _strc:
                    pth, match = thing.check_self(full_part, fix=False)
                    if match:
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
