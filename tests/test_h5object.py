import pytest

from copy import deepcopy

from edges_io.h5 import HDF5RawSpectrum, HDF5StructureExtraKey


def test_extra_key(fastspec_data):
    this = deepcopy(fastspec_data)
    this["meta"]["new_key"] = 2

    with pytest.warns(UserWarning):
        HDF5RawSpectrum.from_data(this)

    HDF5RawSpectrum._require_no_extra = True

    with pytest.raises(HDF5StructureExtraKey):
        HDF5RawSpectrum.from_data(this)

    HDF5RawSpectrum._require_no_extra = False


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
    for key, val in obj.items():
        assert key in obj.keys()
