"""
This module defines the overall file structure and internal contents of the
calibration observations. It does *not* implement any algorithms/methods on that data,
making it easier to separate the algorithms from the data checking/reading.
"""

import glob
import os
import re
import shutil
from abc import ABC, abstractmethod

import h5py
import numpy as np
import read_acq
import toml
from bidict import bidict
from cached_property import cached_property
from scipy import io as sio

from . import utils
from .data import DATA_PATH
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
    def __init__(self, path, fix=False):
        self.path, match = self.check_self(path)

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
    def check_self(path, fix=False):
        pass


class _DataContainer(ABC):
    _content_type = None

    def __init__(self, path, fix=False):
        self.path, match = self.check_self(path, fix)
        self.check_contents(self.path, fix)

        try:
            self._match_dict = match.groupdict()
        except AttributeError:
            try:
                self._match_dict = [m.groupdict() for m in match]
            except TypeError:
                self._match_dict = match
        except Exception:
            raise

        # For a container, if the checks failed, then we ought to bow out now
        if logger.errored:
            logger.errored = False
            raise utils.FileStructureError()

    @classmethod
    @abstractmethod
    def check_self(cls, path, fix=False):
        """Abstract method for checking whether the path is the correct format for
        the DB"""
        pass

    @classmethod
    def check_contents(cls, path, fix=False):
        """Abstract method for checking whether the contents of this container are in
         the correct format for the DB"""
        cls._check_contents_selves(
            path, fix=fix
        )  # Check that everything that *is* there has correct format.
        cls._check_all_files_there(path)  # Check that all necessary files are there.
        # Check that the files that are there have consistent properties, and are also
        # consistent with outside parameters (eg. if year appears on them, they should
        # be consistent with outer years).
        cls._check_file_consistency(path)

    @classmethod
    @abstractmethod
    def _check_all_files_there(cls, path):
        pass

    @classmethod
    @abstractmethod
    def _check_file_consistency(cls, path):
        pass

    @classmethod
    def _check_contents_selves(cls, path, fix=False):
        fls = utils.get_active_files(path)
        for fl in fls:
            if type(cls._content_type) == dict:
                for key, ct in cls._content_type.items():
                    if os.path.basename(fl).startswith(key):
                        content_type = ct
                        break
                else:
                    logger.error(
                        "{} is an extraneous file/folder".format(os.path.basename(fl))
                    )

                    if fix:
                        if os.path.basename(fl) == "Notes.odt":
                            shutil.move(fl, fl.replace("odt", "txt"))
                            fl = fl.replace("odt", "txt")
                            logger.success("Successfully renamed to {}".format(fl))
                            logger.success("Successfully renamed to {}".format(fl))
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
                logger.success("Successfully converted to {}".format(newname))
                shutil.move(os.path.join(root, basename), newpath)
                return newpath, match
            else:
                return None, None

    @classmethod
    def check_self(cls, path, fix=False):
        if type(path) == str:
            path = [path]

        base_fnames = [os.path.basename(fname) for fname in path]
        root = os.path.dirname(os.path.normpath(path[0]))

        matches = []
        for i, basename in enumerate(base_fnames):
            match = re.search(cls.file_pattern, basename)
            if match is None:
                logger.error(
                    "The file {} does not have the correct format for a {}".format(
                        basename, cls.__name__
                    )
                )

                if fix:
                    newname, match = cls._fix(root, basename)
                    path[i] = newname

            if match is not None:
                newname = os.path.basename(path[i])
                groups = match.groupdict()
                if int(groups["run_num"]) < 1:
                    logger.error("The run_num for {} is less than one!".format(newname))
                if not (2010 <= int(groups["year"]) <= 2030):
                    logger.error("The year for {} is a bit strange!".format(newname))
                if not (0 <= int(groups["day"]) <= 366):
                    logger.error(
                        "The day for {} is outside the number of days in a year".format(
                            newname
                        )
                    )
                if not (0 <= int(groups["hour"]) <= 24):
                    logger.error("The hour for {} is outside 0-24!".format(newname))
                if not (0 <= int(groups["minute"]) <= 60):
                    logger.error("The minute for {} is outside 0-60!".format(newname))
                if not (0 <= int(groups["second"]) <= 60):
                    logger.error("The second for {} is outside 0-60!".format(newname))
                if groups["file_format"] not in cls.supported_formats:
                    logger.error(
                        "The file {} is not of a supported format ({}). Got format "
                        "{}".format(
                            newname, cls.supported_formats, groups["file_format"]
                        )
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
                    "No files exist for the load {} with filetype {} on the path: {}".format(
                        load, filetype, direc
                    )
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
                    "No files exist for the load {} for any filetype on that path: {}".format(
                        load, direc
                    )
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

    def read(self):
        """
        Read the files of the object, and concatenate their data.

        Adds the attributes 'p0', 'p1', 'p2' and 'uncalibrated_spectrum'.
        """
        out = {}
        keys = ["p0", "p1", "p2", "Qratio"]
        for fl in self.path:
            this_spec = getattr(self, "_read_" + self.file_format)(fl)

            for key in keys:
                if key not in out:
                    out[key] = this_spec[key]
                else:
                    out[key] = np.concatenate((out[key], this_spec[key]), axis=1)

        return out

    @staticmethod
    def _read_mat(file_name):
        """
        This function loads the antenna temperature and date/time from MAT files.

        Parameters
        ----------
        file_name: str
            The file path to the MAT file.

        Returns
        -------
        2D Uncalibrated Temperature array, or dict of such.
        """
        # loading data and extracting main array
        d = sio.loadmat(file_name)

        # Return dict of all things
        if "Qratio" not in d:
            raise IOError(
                "The file {} is in an old format, and does not have the key Qratio. "
                "Please re-convert it."
            )

        return d

    @staticmethod
    def _read_acq(file_name):
        Q, px = read_acq.decode_file(file_name, progress=False, write_formats=[])
        return {"Qratio": Q.T, "p0": px[0].T, "p1": px[1].T, "p2": px[2].T}

    @staticmethod
    def _read_h5(file_name):
        out = {}
        with h5py.File(file_name, "r") as fl:
            out["Qratio"] = fl["Qratio"][...]
            out["p0"] = fl["p0"][...]
            out["p1"] = fl["p1"][...]
            out["p2"] = fl["p2"][...]

        return out


class Resistance(_SpectrumOrResistance):
    """
    An object representing a resistance measurement (and its structure).
    """

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
            return self._data
        except AttributeError:
            data = np.genfromtxt(
                self.path,
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

            if self.store_data:
                self._data = data

            return data

    @property
    def raw_data(self):
        """The raw csv data of the Resistance measurement.

        Note that this is only cached in memory if `store_data` is True, otherwise
        the data is re-read from disk each time `raw_data` is accessed.
        """
        return self.read()


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
    def _check_all_files_there(cls, path):
        # Just need to check for the loads.
        for name, load in LOAD_ALIASES.items():
            if not glob.glob(os.path.join(path, load + "_*")):
                logger.error(
                    "{} does not contain any files for load {}".format(
                        cls.__name__, load
                    )
                )

    @classmethod
    def get_all_load_names(cls, path):
        """Get all load names found in the Spectra directory"""
        fls = utils.get_active_files(path)
        return {os.path.basename(fl).split("_")[0] for fl in fls}

    @classmethod
    def get_simulator_names(cls, path):
        load_names = cls.get_all_load_names(path)
        return {name for name in load_names if name in ANTENNA_SIMULATORS}

    @classmethod
    def _check_file_consistency(cls, path):
        fls = utils.get_active_files(path)

        groups = [
            re.search(cls._content_type.file_pattern, fl).groupdict() for fl in fls
        ]

        # Ensure all years are the same
        for fl, group in zip(fls, groups):
            if group["year"] != groups[0]["year"]:
                logger.error(
                    "All years must be the same in a Spectra folder, but {} "
                    "was not".format(fl)
                )

        # Ensure days are close-ish
        days = [int(group["day"]) for group in groups]
        if max(days) - min(days) > 30:
            logger.error(
                "Observation days are suspiciously far apart for {}".format(path)
            )

    def read_all(self):
        """Read all spectra"""
        return {name: getattr(self, name).read() for name in LOAD_ALIASES}


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

        #  loading data
        d = np.genfromtxt(path_filename, skip_header=comment_rows)

        return d, flag


class _S11SubDir(_DataContainer):
    STANDARD_NAMES = S1P.POSSIBLE_KINDS
    _content_type = S1P
    folder_pattern = None

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
    def _check_all_files_there(cls, path):
        for name in cls.STANDARD_NAMES:
            if not glob.glob(os.path.join(path, name + "??.s1p")):
                logger.error("No {} standard found in {}".format(name, path))

    @classmethod
    def _check_file_consistency(cls, path):
        pass

    def _get_max_run_num(self):
        return max(
            int(re.match(S1P.file_pattern, os.path.basename(fl)).group("run_num"))
            for fl in self.active_contents
        )


class LoadS11(_S11SubDir):
    STANDARD_NAMES = ["Open", "Short", "Match", "External"]
    folder_pattern = "(?P<load_name>{})$".format("|".join(LOAD_ALIASES.values()))

    def __init__(self, direc, run_num=None, fix=False):
        super().__init__(direc, run_num, fix)
        self.load_name = LOAD_ALIASES.inverse.get(
            self._match_dict["load_name"], self._match_dict["load_name"]
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

    def __init__(self, path, repeat_num=None, run_num=None, fix=False):
        """Class representing the entire S11 subdirectory of an observation

        Parameters
        ----------
        path : str
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
        fls = utils.get_active_files(os.path.join(path))
        fls = [fl for fl in fls if kind in fl]
        rep_nums = [int(fl[-2:]) for fl in fls]
        return max(rep_nums)

    @classmethod
    def check_self(cls, path, fix=False):
        logger.structure("Checking S11 folder contents at {}".format(path))

        if not os.path.exists(path):
            logger.error("This path does not exist: {}".format(path))

        if os.path.basename(path) != "S11":
            logger.error("The S11 folder should be called S11")

        return path, True

    @classmethod
    def _check_all_files_there(cls, path):
        for load in LOAD_ALIASES.values():
            if not glob.glob(os.path.join(path, load)):
                logger.error("No {} S11 directory found!".format(load))

        for other in ["SwitchingState", "ReceiverReading"]:
            if not glob.glob(os.path.join(path, other + "??")):
                logger.error("No {} S11 directory found!".format(other))

    @classmethod
    def _check_file_consistency(cls, path):
        simulators = cls.get_simulator_names(path)
        if simulators:
            logger.info(
                "Found the following Antenna Simulators in S11: {}".format(
                    ",".join(simulators)
                )
            )
        else:
            logger.info("No Antenna Simulators in S11.")

    @classmethod
    def get_simulator_names(cls, path):
        fls = utils.get_active_files(path)
        return {
            os.path.basename(fl)
            for fl in fls
            if any(os.path.basename(fl).startswith(k) for k in ANTENNA_SIMULATORS)
        }


class CalibrationObservation(_DataContainer):
    file_pattern = (
        r"^Receiver(?P<rcv_num>\d{2})_(?P<year>\d{4})_(?P<month>\d{2})_(?P<day>\d{2})_("
        r"?P<freq_low>\d{3})_to_(?P<freq_hi>\d{3})MHz/(?P<temp>\d{2})C$"
    )
    _content_type = {"S11": S11Dir, "Spectra": Spectra, "Resistance": Resistances}

    def __init__(self, path, ambient_temp=25, run_num=None, repeat_num=None, fix=False):
        """
        A class defining a full calibration observation, with all Spectra, Resistance
        and S11 files necessary to do a single analysis.
        """
        if ambient_temp not in [15, 25, 35]:
            raise ValueError("ambient temp must be one of 15, 25, 35!")

        path = os.path.join(os.path.abspath(path), "{}C".format(ambient_temp))

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
            os.path.join(self.path, "Spectra"),
            run_num=run_nums.get("Spectra", None),
            fix=fix,
        )
        self.resistance = Resistances(
            os.path.join(self.path, "Resistance"),
            run_num=run_nums.get("Resistance", None),
            fix=fix,
        )
        self.s11 = S11Dir(
            os.path.join(self.path, "S11"),
            run_num=run_nums.get("S11", None),
            repeat_num=repeat_num,
            fix=fix,
        )

        self.simulator_names = self.get_simulator_names(self.path)

    @classmethod
    def check_self(cls, path, fix=False):
        logger.structure("Checking root folder: {}".format(path))

        if not os.path.exists(path):
            raise IOError("The path {} does not exist!".format(path))

        path = os.path.normpath(path)
        base = os.path.dirname(os.path.dirname(path))
        name = os.path.join(
            os.path.basename(os.path.dirname(path)), os.path.basename(path)
        )

        match = re.search(cls.file_pattern, name)

        if match is None:
            logger.error(
                "Calibration Observation directory name is in the wrong format!"
            )

            if fix:

                bad_pattern = (
                    r"^Receiver(\d{1,2})_(\d{4})_(\d{1,2})_(\d{1,2})_(\d{2,3})_to_(\d{"
                    r"2,3})MHz/(?P<temp>\d{2})C$"
                )

                match = re.search(bad_pattern, name)

                if match is None:
                    bad_pattern = (
                        r"Receiver(?P<rcv_num>\d{2})_(?P<year>\d{4})_(?P<month>\d{2})_(?P<day>\d{2})_("
                        r"?P<freq_low>\d{3})_to_(?P<freq_hi>\d{3})_MHz/(?P<temp>\d{2})C$"
                    )
                    match = re.search(bad_pattern, name)

                if match is not None:
                    newname = "Receiver{:0>2}_{}_{:0>2}_{:0>2}_{:0>3}_to_{:0>3}MHz/{}C".format(
                        *match.groups()
                    )
                    shutil.move(
                        os.path.normpath(os.path.dirname(path)),
                        os.path.join(base, os.path.dirname(newname)),
                    )

                    # # If top-level directory is now empty, remove it.
                    # if not glob.glob(os.path.join(os.path.dirname(os.path.normpath(path))), "*"):
                    #     os.rmdir(os.path.dirname(os.path.normpath(path)))

                    name = newname
                    path = os.path.join(base, name)
                    logger.success("Successfully renamed to {}".format(newname))
                else:
                    logger.warning("Failed to fix the name scheme")

        if match is not None:
            groups = match.groupdict()
            if int(groups["rcv_num"]) < 1:
                logger.error("Unknown receiver number: {}".format(groups["rcv_num"]))
            if not (2010 <= int(groups["year"]) <= 2030):
                logger.error("Unknown year: {}".format(groups["year"]))
            if not (1 <= int(groups["month"]) <= 12):
                logger.error("Unknown month: {}".format(groups["month"]))
            if not (1 <= int(groups["day"]) <= 31):
                logger.error("Unknown day: {}".format(groups["day"]))
            if not (1 <= int(groups["freq_low"]) <= 300):
                logger.error("Low frequency is weird: {}".format(groups["freq_low"]))
            if not (1 <= int(groups["freq_hi"]) <= 300):
                logger.error("High frequency is weird: {}".format(groups["high_hi"]))
            if not int(groups["freq_low"]) < int(groups["freq_hi"]):
                logger.error(
                    "Low frequency > High Frequency: {} > {}".format(
                        groups["freq_low"]
                    ),
                    groups["freq_hi"],
                )

            logger.info("Calibration Observation Metadata: {}".format(groups))

        logger.debug("\tReturning path={}".format(path))

        return path, match

    @classmethod
    def _check_all_files_there(cls, path):
        for folder in ["S11", "Spectra", "Resistance"]:
            if not os.path.exists(os.path.join(path, folder)):
                logger.error("No {} folder in observation!".format(folder))

    @classmethod
    def _check_file_consistency(cls, path):
        cls.get_simulator_names(
            path
        )  # checks whether simulators are the same between each.

    @classmethod
    def get_simulator_names(cls, path):
        # Go through the subdirectories and check their simulators
        dct = {
            name: tuple(sorted(kls.get_simulator_names(os.path.join(path, name))))
            for name, kls in cls._content_type.items()
        }

        # If any list of simulators is not the same as the others, make an error.
        if len(set(dct.values())) != 1:
            logger.error(
                "Antenna Simulators do not match in all subdirectories. Got {}".format(
                    dct
                )
            )

        return set(list(dct.values())[0])

    def read_all(self):
        """Read all spectra and resistance files."""
        self.spectra.read_all()
        self.resistance.read_all()
