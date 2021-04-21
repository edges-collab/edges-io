import attr
import contextlib
import h5py
import numpy as np
import warnings
from abc import ABCMeta, abstractmethod
from datetime import datetime
from pathlib import Path

from . import __version__, utils


class HDF5StructureError(Exception):
    pass


class HDF5StructureValidationError(HDF5StructureError):
    pass


class HDF5StructureExtraKey(HDF5StructureError):
    pass


@attr.s
class _HDF5Part(metaclass=ABCMeta):
    filename = attr.ib(default=None, converter=lambda x: x if x is None else Path(x))
    group_path = attr.ib(default="", converter=str)

    def __attrs_post_init__(self):
        self.__memcache__ = {}
        self.__fl_inst = None

    @contextlib.contextmanager
    def open(self) -> h5py.Group:
        """Context manager for opening up the file.

        Yields
        ------
        grp : :class:`h5py.Group`
            The h5py Group corresponding to this instance.
        """
        grp = self._fl_instance

        if self.group_path:
            for bit in self.group_path.split("."):
                grp = grp[bit]

        yield grp

    @property
    def _fl_instance(self):
        if not self.filename:
            raise IOError(
                "This object has no associated file. You can define one with the write() method."
            )

        if self.__fl_inst is None:
            self.__fl_inst = h5py.File(self.filename, "r")

        return self.__fl_inst

    def __getstate__(self):
        """Prepare class for pickling. HDF5 files are not pickleable!"""
        return {
            key: val
            for key, val in self.__dict__.items()
            if not key.endswith("__fl_inst")
        }

    def __contains__(self, item):
        return item in list(self.keys())

    def __getitem__(self, item):
        if item in self.__memcache__:
            return self.__memcache__[item]

        with self.open() as fl:
            if item in ("attrs", "meta"):
                out = dict(fl.attrs)
                for k, v in out.items():
                    if isinstance(v, str) and v == "none":
                        out[k] = None
            elif item not in fl:
                raise KeyError(
                    f"'{item}' is not a valid part of {self.__class__.__name__}."
                    f" Valid keys: {self.keys()}"
                )
            elif isinstance(fl[item], h5py.Group):
                if not isinstance(self._structure[item], dict):
                    raise HDF5StructureValidationError(
                        f"item {item} has structure {self._structure[item]}, but must be dict."
                    )

                gp = self.group_path + "." if self.group_path else ""
                out = _HDF5Group(
                    filename=self.filename,
                    structure=self._structure[item],
                    group_path=gp + item,
                )

            elif isinstance(fl[item], h5py.Dataset):
                out = fl[item][...]
            else:
                raise NotImplementedError("that item is not supported yet.")

        # Save the key to the cache.
        self.__memcache__[item] = out

        return out

    def keys(self):
        with self.open() as fl:
            yield from fl.keys()

    def items(self):
        for k in self.keys():
            yield k, self[k]

    def clear(self, keys=None):
        """Clear all items from memory loaded from file."""
        if keys is None:
            self.__memcache__ = {}
        else:
            for key in keys:
                del self.__memcache__[key]

    def cached_keys(self):
        """All of the keys that have already been read into cache."""
        return self.__memcache__.keys()

    @property
    def meta(self):
        """Metadata of the object."""
        return self["meta"]


@attr.s
class HDF5Object(_HDF5Part):
    """
    An object that provides a transparent wrapper of a HDF5 file.

    Creation of this object can be done in two ways: either by passing a filename
    to wrap, or by using ``.from_data``, in which case you must pass all the data to it,
    which it can write in the correct format.

    This class exists to be subclassed. Subclasses should define the attribute
    ``_structure``, which defines the layout of the underlying file. The attribute
    ``_require_all`` sets whether checks on the file will fail if not all keys in the
    structure are present in the file. Conversely ``_require_no_extra`` sets whether
    it will fail if extra keys are present.

    Parameters
    ----------
    filename : str or Path
        The filename of the HDF5 file to wrap.
    require_all : bool, optional
        Over-ride the class attribute requiring all the structure to exist in file.
    require_no_extra : bool, optional
        Over-ride the class attribute requiring no extra data to exist in file.

    Notes
    -----
    Accessing data is very similar to just using the `h5py.File`,
    except that the object is able to check the structure of the file, and is slightly
    more convenient.

    Note that an `.open` method exists which returns an open ``h5py.File`` object. This
    is a context manager so you can do things like::

        with obj.open() as fl:
            val = fl.attrs['key']
    """

    _require_no_extra = False
    default_root = Path(".")
    _structure = None

    require_no_extra = attr.ib(default=_require_no_extra, converter=bool)

    def __attrs_post_init__(self):
        super().__attrs_post_init__()

        if self.filename and self.filename.exists():
            self.check(self.filename, self.require_no_extra)

    def __init_subclass__(cls, **kwargs):
        super().__init_subclass__()
        if "meta" not in cls._structure:
            cls._structure["meta"] = {}

        for k, v in cls._get_extra_meta().items():
            if k not in cls._structure["meta"]:
                cls._structure["meta"][k] = None

    @classmethod
    def from_data(cls, data, validate=True, **kwargs):
        inst = cls(**kwargs)

        if "meta" not in data:
            data["meta"] = {}

        data["meta"].update(cls._get_extra_meta())

        if validate:
            false_if_extra = kwargs.get("require_no_extra", cls._require_no_extra)

            try:
                cls._checkgrp(data, cls._structure)
            except HDF5StructureExtraKey as e:
                if false_if_extra:
                    raise HDF5StructureExtraKey(
                        f"Data had extra key(s)! Extras: {str(e).split(':')[-1]}"
                    )
                else:
                    warnings.warn(
                        f"Data had extra key! Extras: {str(e).split(':')[-1]}"
                    )

        inst.__memcache__ = data
        return inst

    @classmethod
    def _get_extra_meta(cls):
        return {
            "write_time": datetime.now().strftime(
                datetime.now().strftime("%Y-%M-%D:%H:%M:%S")
            ),
            "edges_io_version": __version__,
            "object_name": cls.__name__,
        }

    def write(self, filename=None, clobber=False):
        if filename is None and self.filename is None:
            raise ValueError(
                "You need to pass a filename since there is no instance filename."
            )

        filename = Path(filename or self.filename)

        if not filename.is_absolute():
            filename = self.default_root / filename

        if self.filename is None:
            self.filename = filename

        if filename.exists() and not clobber:
            raise FileExistsError(f"file {filename} already exists!")

        def _write(grp, struct, cache):
            for k, v in cache.items():
                try:
                    if isinstance(v, dict):
                        g = grp.attrs if k in ["meta", "attrs"] else grp.create_group(k)
                        _write(g, struct[k], v)
                    else:
                        if v is None:
                            v = "none"
                        elif isinstance(v, Path):
                            v = str(v)

                        grp[k] = v

                except TypeError:
                    raise TypeError(
                        f"For key '{k}' in class '{self.__class__.__name__}', type '"
                        f"{type(v)}' is not allowed in HDF5."
                    )

        to_write = self.__memcache__

        if "meta" not in to_write:
            to_write["meta"] = {}

        to_write["meta"].update(self._get_extra_meta())

        if not filename.parent.exists():
            filename.parent.mkdir(parents=True)

        with h5py.File(filename, "w") as fl:
            _write(fl, self._structure, to_write)

    @classmethod
    def _checkgrp(cls, grp, strc):
        for k, v in strc.items():
            # We treat 'meta' as synonymous with 'attrs'
            if k == "meta" and k not in grp:
                k = "attrs"

            if (
                k not in grp
                and k != "attrs"
                and v != "optional"
                and not getattr(v, "optional", False)
            ):
                raise TypeError(f"Non-optional key '{k}' not in {grp}")
            elif k == "attrs":
                if isinstance(grp, (h5py.Group, h5py.File)):
                    cls._checkgrp(grp.attrs, v)
                else:
                    cls._checkgrp(grp[k], v)
            elif isinstance(v, dict):
                cls._checkgrp(grp[k], v)
            elif not (v is None or v == "optional" or v(grp.get(k, None))):
                raise HDF5StructureValidationError(
                    f"key {k} in {grp} failed its validation. Type: {type(grp[k])}"
                )

        # Ensure there's no extra keys in the group
        if len(strc) < len(grp.keys()):
            extras = [k for k in grp.keys() if k not in strc]
            raise HDF5StructureExtraKey(f"Extra keys found in {grp}: {extras}")

    @classmethod
    def check(cls, filename, false_if_extra=None):
        false_if_extra = false_if_extra or cls._require_no_extra

        if not cls._structure:
            return True

        with h5py.File(filename, "r") as fl:
            try:
                cls._checkgrp(fl, cls._structure)
            except HDF5StructureExtraKey as e:
                if false_if_extra:
                    raise e
                else:
                    warnings.warn(f"{e}. Filename={filename}. ")


@attr.s
class _HDF5Group(_HDF5Part):
    """Similar to HDF5Object, but pointing to a Group within it."""

    _structure = attr.ib(factory=dict, converter=dict)


class HDF5RawSpectrum(HDF5Object):
    _require_no_extra = False

    _structure = {
        "meta": {
            "fastspec_version": utils.optional(utils.isstringish),
            "start": utils.optional(utils.isintish),
            "stop": utils.optional(utils.isintish),
            "site": utils.optional(utils.isstringish),
            "instrument": utils.optional(utils.isstringish),
            "switch_delay": utils.optional(utils.isfloatish),
            "input_channel": utils.optional(utils.isintish),
            "voltage_range": utils.optional(utils.isnumeric),
            "samples_per_accumulation": utils.optional(utils.isintish),
            "acquisition_rate": utils.optional(utils.isnumeric),
            "num_channels": utils.optional(utils.isintish),
            "num_taps": utils.optional(utils.isintish),
            "window_function_id": utils.optional(utils.isintish),
            "num_fft_threads": utils.optional(utils.isintish),
            "num_fft_buffers": utils.optional(utils.isintish),
            "stop_cycles": utils.optional(utils.isintish),
            "stop_seconds": utils.optional(utils.isfloatish),
            "stop_time": "optional",
            "edges_io_version": utils.isstringish,
            "object_name": utils.isstringish,
            "write_time": utils.isstringish,
            "show": "optional",  # From here down, we just include them as optional for backwards compat.
            "hide": "optional",
            "kill": "optional",
            "help": "optional",
            "inifile": "optional",
            "datadir": "optional",
            "output_file": "optional",
            "switch_io_port": "optional",
            "samples_per_transfer": "optional",
            "show_plots": "optional",
            "plot_bin": "optional",
            "resolution": "optional",
            "temperature": "optional",
            "nblk": "optional",
            "nfreq": "optional",
            "freq_min": "optional",
            "freq_max": "optional",
            "freq_res": "optional",
            "n_file_lines": "optional",
        },
        "spectra": {
            "p0": lambda x: (x.ndim == 2 and x.dtype == float),
            "p1": lambda x: (x.ndim == 2 and x.dtype == float),
            "p2": lambda x: (x.ndim == 2 and x.dtype == float),
            "Q": lambda x: (x.ndim == 2 and x.dtype == float),
        },
        "freq_ancillary": {"frequencies": lambda x: (x.ndim == 1 and x.dtype == float)},
        "time_ancillary": {
            "times": lambda x: (x.ndim == 1 and x.dtype == "|S17"),
            "adcmax": lambda x: (
                x.ndim == 2 and x.shape[1] == 3 and x.dtype in (float, np.float32)
            ),
            "adcmin": lambda x: (
                x.ndim == 2 and x.shape[1] == 3 and x.dtype in (float, np.float32)
            ),
            "data_drops": "optional",
        },
    }
