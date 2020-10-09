import attr
import contextlib
import h5py
import numpy as np
import warnings
from datetime import datetime
from pathlib import Path

from . import __version__


class HDF5StructureError(Exception):
    pass


class HDF5StructureValidationError(HDF5StructureError):
    pass


class HDF5StructureExtraKey(HDF5StructureError):
    pass


@attr.s
class HDF5Object:
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
    lazy : bool, optional
        Whether to be lazy in loading the data in the file into the structure.

    Notes
    -----
    With ``lazy=True``, accessing data is very similar to just using the `h5py.File`,
    except that the object is able to check the structure of the file, and is slightly
    more convenient.

    Note that an `.open` method exists which returns an open ``h5py.File`` object. This
    is a context manager so you can do things like::

        with obj.open() as fl:
            val = fl.attrs['key']

    """

    _structure = None
    _require_no_extra = False
    default_root = Path(".")

    filename = attr.ib(default=None, converter=lambda x: x if x is None else Path(x))
    require_no_extra = attr.ib(default=_require_no_extra, converter=bool)
    lazy = attr.ib(default=True, converter=bool)

    def __attrs_post_init__(self):
        self.__memcache__ = {}

        if self.filename and self.filename.exists():
            self.check(self.filename, self.require_no_extra)

        if not self.lazy:
            self.load_all(self.filename)

    @classmethod
    def from_data(cls, data, **kwargs):
        inst = cls(**kwargs)

        false_if_extra = kwargs.get("require_no_extra", cls._require_no_extra)

        try:
            cls._checkgrp(data, cls._structure)
        except HDF5StructureExtraKey as e:
            if false_if_extra:
                raise HDF5StructureExtraKey(
                    f"Data had extra key(s)! Extras: {str(e).split(':')[-1]}"
                )
            else:
                warnings.warn(f"Data had extra key! Extras: {str(e).split(':')[-1]}")

        inst.__memcache__ = data
        return inst

    @contextlib.contextmanager
    def open(self, mode: str = "r") -> h5py.File:
        """Context manager to open the file.

        Parameters
        ----------
        mode : str
            The read/write mode to open the file in.

        Yields
        ------
        fl : `h5py.File` instance
            An instance of the open file.
        """
        fl = h5py.File(self.filename, mode=mode)
        yield fl
        fl.close()

    def load(self, key: str) -> [dict, h5py.Dataset, h5py.Group]:
        """Load key from file into memory and keep it cached in memory.

        Parameters
        ----------
        key : str
            The key to load in.

        Returns
        -------
        out : dict, :class:`h5py.Group` or :class:`h5py.Dataset`
            The dataset or group of such datasets that were loaded.

        Notes
        -----
        If you would like to load in a key further down the
        heirarchy *without* loading in everything above it (i.e. just getting that
        particular bit of memory mapped in), then use item-getting to get to the
        parent of this key. This will load *everything* into memory below the given key.

        Examples
        --------
        Consider the object

        >>> h5obj = HDF5Object("filename.h5")

        which has a top-level dataset called "top-level". To load in this dataset to
        memory:

        >>> top_level = h5obj.load('top-level')

        You can use the ``top_level`` array, or access the same memory in a fast way
        via

        >>> h5obj['top-level']  # returns array from memory (rather than loading from file)

        Now, consider a group at the top level called "group", under which is a dataset
        called ``dataset`` and another dataset called ``superfluous`` (which we don't
        need right now). To load ``dataset`` into memory without loading ``superfluous``,
        use

        >>> dataset = h5obj['group'].load('dataset')

        Alternatively, if you just want to load the whole group in and have access to
        both datasets from memory:

        >>> group = h5obj.load('group')
        >>> dataset = group['dataset']
        """
        # Load up this key into memory.
        # This checks the memcache first.
        out = self[key]
        try:
            # If the key was a group, it hasn't loaded any data yet, so load that.
            out.load_all()
        except AttributeError:
            pass

        # Save the key to the cache.
        if key not in self.__memcache__:
            self.__memcache__[key] = out
        return out

    def load_all(self, filename=None):
        if filename and not self.filename:
            self.filename = filename

        filename = filename or self.filename

        if not filename:
            raise ValueError("You need to provide a filename to load!")

        for k in self._structure.keys():
            self.load(k)

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

            if k not in grp and k != "attrs" and v != "optional":
                raise TypeError(f"Non-optional key '{k}' not in {grp}")
            elif k == "attrs":
                if isinstance(grp, (h5py.Group, h5py.File)):
                    cls._checkgrp(grp.attrs, v)
                else:
                    cls._checkgrp(grp[k], v)
            elif isinstance(v, dict):
                cls._checkgrp(grp[k], v)
            elif not (v is None or v == "optional" or v(grp[k])):
                raise HDF5StructureValidationError(
                    f"key {k} in {grp} failed its validation. Type: {type(grp[k])}"
                )

        # Ensure there's no extra keys in the group
        if len(strc) < len(grp.keys()):
            extras = [k for k in grp.keys() if k not in strc]
            raise HDF5StructureExtraKey(f"Extra keys found in the file: {extras}")

    @classmethod
    def check(cls, filename, false_if_extra=None):
        false_if_extra = false_if_extra or cls._require_no_extra

        if not cls._structure:
            return True

        with h5py.File(filename, "r") as fl:
            try:
                cls._checkgrp(fl, cls._structure)
            except HDF5StructureExtraKey:
                if false_if_extra:
                    raise HDF5StructureExtraKey(f"Extra key found in {filename}")
                else:
                    warnings.warn(f"Extra key found in {filename}")

    def __contains__(self, item):
        return item in self.keys()

    def __getitem__(self, item):
        if item in self.__memcache__:
            return self.__memcache__[item]

        if item not in self._structure:
            raise KeyError(
                f"'{item}' is not a valid part of {self.__class__.__name__}."
                f" Valid keys: {self.keys()}"
            )

        with h5py.File(self.filename, "r") as fl:
            if item in ("attrs", "meta"):
                out = dict(fl.attrs)
                for k, v in out.items():
                    if isinstance(v, str) and v == "none":
                        out[k] = None
            elif isinstance(fl[item], h5py.Group):
                out = _HDF5Group(self.filename, self._structure[item], item)

                # If it's a group, we *do* want to add it to cache, since otherwise
                # its own cache is lost.
                self.__memcache__[item] = out
            elif isinstance(fl[item], h5py.Dataset):
                out = fl[item][...]
            else:
                raise NotImplementedError("that item is not supported yet.")

        return out

    def keys(self):
        return self._structure.keys()

    def items(self):
        for k in self.keys():
            yield k, self[k]


@attr.s
class _HDF5Group:
    """Similar to HDF5Object, but pointing to a Group within it."""

    filename = attr.ib(converter=Path, validator=lambda x, att, val: val.exists())
    structure = attr.ib(converter=dict)
    group_path = attr.ib(converter=str)
    lazy = attr.ib(default=True, converter=bool)

    def __attrs_post_init__(self):
        self.__memcache__ = {}

        if not self.lazy:
            self.load_all()

    def load(self, key: str):
        return HDF5Object.load(self, key)

    def load_all(self):
        """Read the whole file into memory at once."""
        for k in self.structure:
            self.load(k)

    @contextlib.contextmanager
    def open(self) -> h5py.Group:
        """Context manager for opening up the file and getting this group.

        Yields
        ------
        grp : :class:`h5py.Group`
            The h5py Group corresponding to this instance.
        """
        fl = h5py.File(self.filename, "r")
        grp = fl

        for bit in self.group_path.split("."):
            grp = grp[bit]

        yield grp

        fl.close()

    def __getitem__(self, item):
        if item in self.__memcache__:
            return self.__memcache__[item]

        with self.open() as fl:
            if item in ("attrs", "meta"):
                out = dict(fl.attrs)
            elif isinstance(fl[item], h5py.Group):
                out = _HDF5Group(self.filename, item)
            elif isinstance(fl[item], h5py.Dataset):
                out = fl[item][...]
            else:
                raise NotImplementedError("that item is not supported yet.")

        return out

    def keys(self):
        return self.structure.keys()

    def items(self):
        for k in self.keys():
            yield k, self[k]


class HDF5RawSpectrum(HDF5Object):
    _require_no_extra = False

    _structure = {
        "meta": {
            "fastspec_version": "optional",
            "start": "optional",  # lambda x: isinstance(x, (int, np.int64)),
            "stop": "optional",  # lambda x: isinstance(x, (int, np.int64)),
            "site": "optional",  # lambda x: isinstance(x, str),
            "instrument": "optional",  # lambda x: isinstance(x, str),
            "switch_io_port": "optional",  # lambda x: isinstance(x, (int, np.int64)),
            "switch_delay": "optional",  # lambda x: isinstance(x, np.float),
            "input_channel": "optional",  # lambda x: isinstance(x, (int, np.int64)),
            "voltage_range": "optional",  # lambda x: isinstance(x, (int, np.int64)),
            "samples_per_accumulation": "optional",  # lambda x: isinstance(x, (int, np.int64)),
            "acquisition_rate": "optional",  # lambda x: isinstance(x, (int, np.int64)),
            "num_channels": "optional",  # lambda x: isinstance(x, (int, np.int64)),
            "num_taps": "optional",  # lambda x: isinstance(x, (int, np.int64)),
            "window_function_id": "optional",  # lambda x: isinstance(x, (int, np.int64)),
            "num_fft_threads": "optional",  # lambda x: isinstance(x, (int, np.int64)),
            "num_fft_buffers": "optional",  # lambda x: isinstance(x, (int, np.int64)),
            "stop_cycles": "optional",  # lambda x: isinstance(x, (int, np.int64)),
            "stop_seconds": "optional",  # lambda x: isinstance(x, np.float),
            "stop_time": "optional",  # lambda x: isinstance(x, str),
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
