import pytest

import dill as pickle
import h5py
import inspect
import numpy as np
from copy import deepcopy
from pathlib import Path

from edges_io.h5 import (
    HDF5Object,
    HDF5RawSpectrum,
    HDF5StructureExtraKey,
    HDF5StructureValidationError,
)


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


def test_write_no_fname(fastspec_data):
    obj = HDF5RawSpectrum.from_data(fastspec_data)

    with pytest.raises(ValueError):
        obj.write()


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
    obj["spectra"]["p0"]

    assert "p0" in obj.__memcache__["spectra"].__memcache__
    obj["spectra"].clear()
    assert "p0" not in obj.__memcache__["spectra"].__memcache__

    obj.clear(["spectra"])
    assert "spectra" not in obj.__memcache__


def test_getitem(fastspec_spectrum_fl):
    obj = HDF5RawSpectrum(fastspec_spectrum_fl)
    assert "spectra" in obj
    assert obj["spectra"]["p0"].ndim == 2
    assert isinstance(obj["meta"], dict)
    assert isinstance(obj["attrs"], dict)
    assert obj["attrs"] == obj["meta"]

    with pytest.raises(KeyError):
        obj["not_existent"]


def test_pickling(fastspec_spectrum_fl):
    obj = HDF5RawSpectrum(fastspec_spectrum_fl)
    # ensure we load up the file instance
    obj._fl_instance

    # see if we can still pickle it...
    pickle.dumps(obj)


def test_bad_existing_h5(tmpdir: Path):
    class Bad(HDF5Object):
        _structure = {"data": lambda x: isinstance(x, np.ndarray)}

    with h5py.File(tmpdir / "bad.h5", "w") as fl:
        grp = fl.create_group("data")
        grp.attrs["bad_key"] = True
        grp["data"] = np.linspace(0, 1, 10)

    with pytest.raises(HDF5StructureValidationError):
        Bad(tmpdir / "bad.h5")


def test_h5_hierarchical(tmpdir: Path):
    class Example(HDF5Object):
        _structure = {
            "this": {"that": {"the_other": {"key": lambda x: x.shape == (10,)}}}
        }

    ex = Example.from_data({"this": {"that": {"the_other": {"key": np.zeros(10)}}}})

    assert ex["this"]["that"]["the_other"]["key"].shape == (10,)

    ex.write(tmpdir / "tmp_hierarchical.h5")

    ex2 = Example(tmpdir / "tmp_hierarchical.h5")

    assert ex2["this"]["that"]["the_other"]["key"].shape == (10,)

    assert isinstance(ex2.meta, dict)
