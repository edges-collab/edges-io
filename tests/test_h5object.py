import pytest

import h5py
import inspect
import numpy as np
from copy import deepcopy
from pathlib import Path

from edges_io.h5 import HDF5Object, HDF5RawSpectrum, HDF5StructureExtraKey


def test_extra_key(fastspec_data, fastspec_spectrum_fl):
    this = deepcopy(fastspec_data)
    this["meta"]["new_key"] = 2
    this["spectra"]["new_spectrum"] = np.zeros(20)

    with pytest.warns(UserWarning):
        obj = HDF5RawSpectrum.from_data(this)

    # Ensure we can still access 'new_spectrum' even though it's not in the _structure
    assert len(obj["spectra"]["new_spectrum"]) == 20

    # But we can't get stuff that doesn't exist at all.
    with pytest.raises(KeyError):
        obj["spectra"]["non-existent"]

    with pytest.raises(IOError):
        obj["non-existent"]

    HDF5RawSpectrum._require_no_extra = True

    with pytest.raises(HDF5StructureExtraKey):
        HDF5RawSpectrum.from_data(this, require_no_extra=True)

    HDF5RawSpectrum._require_no_extra = False

    this = HDF5RawSpectrum(fastspec_spectrum_fl)
    with pytest.raises(KeyError):
        this["non-existent"]


def test_h5_open(fastspec_spectrum_fl):
    obj = HDF5RawSpectrum(fastspec_spectrum_fl)
    with obj.open():
        assert "spectra" in obj.keys()


def test_load_all_no_fname(fastspec_data):
    obj = HDF5RawSpectrum.from_data(fastspec_data)
    with pytest.raises(ValueError):
        obj.load_all()

    with pytest.raises(ValueError):
        obj.write()


def test_load_all(fastspec_spectrum_fl):
    obj = HDF5RawSpectrum(fastspec_spectrum_fl)
    obj.load_all()
    assert "spectra" in obj.__memcache__


def test_access_nonexistent(fastspec_spectrum_fl):
    obj = HDF5RawSpectrum(fastspec_spectrum_fl)
    with pytest.raises(KeyError):
        obj["nonexistent"]


def test_access_attrs(fastspec_spectrum_fl):
    obj = HDF5RawSpectrum(fastspec_spectrum_fl)
    assert "stop" in obj["meta"]


def test_get_items(fastspec_spectrum_fl):
    obj = HDF5RawSpectrum(fastspec_spectrum_fl)
    assert inspect.isgeneratorfunction(obj.items)
    for key, val in obj.items():
        assert key in obj.keys()


def test_read_none(tmpdir: Path):
    fname = tmpdir / "tmp_file.h5"

    with h5py.File(fname, "w") as fl:
        fl.attrs["key"] = "none"

    obj = HDF5Object(fname)
    assert obj["meta"]["key"] is None
    assert obj["attrs"]["key"] is None


def test_read_group_meta(fastspec_spectrum_fl):
    obj = HDF5RawSpectrum(fastspec_spectrum_fl)
    assert obj["spectra"]["meta"] == {}

    with pytest.raises(KeyError):
        obj["spectra"]["non-existent"]

    assert len(list(obj["spectra"].keys())) == 4
    for key, val in obj["spectra"].items():
        assert isinstance(val, np.ndarray)


def test_clear(fastspec_spectrum_fl):
    obj = HDF5RawSpectrum(fastspec_spectrum_fl)
    obj.load_all()
    assert "p0" in obj.__memcache__["spectra"].__memcache__
    obj["spectra"].clear()
    assert "p0" not in obj.__memcache__["spectra"].__memcache__
    obj.clear()
    assert "spectra" not in obj.__memcache__
