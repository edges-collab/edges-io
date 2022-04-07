import pytest

import numpy as np
from pathlib import Path

from edges_io.h5 import HDF5RawSpectrum
from edges_io.io import Spectrum


def check_obj(obj):
    assert "spectra" in obj
    assert "freq_ancillary" in obj
    assert "time_ancillary" in obj
    assert "start" in obj.meta

    # Look at an item
    b = obj["spectra"]["Q"]
    assert b.ndim == 2
    assert "Q" in obj["spectra"].__memcache__
    assert "p0" not in obj["spectra"].__memcache__

    t = obj["time_ancillary"]["times"]
    if t.ndim == 2:
        tt = t[:, 0]
    else:
        tt = t
    assert obj.convert_times(tt) == obj.get_times()

    if t.ndim == 1:
        with pytest.warns(UserWarning, match="Cannot read times for swpos > 0"):
            obj.get_times(swpos=1)


def test_hdf5rawspectrum(fastspec_spectrum_fl):
    obj = HDF5RawSpectrum(fastspec_spectrum_fl)
    check_obj(obj)


def test_hdf5rawspectrum_2dim_time(fastspec_spectrum_fl_2dim_time):
    obj = HDF5RawSpectrum(fastspec_spectrum_fl_2dim_time)
    check_obj(obj)


def test_io_read(fastspec_spectrum_fl):
    spec = Spectrum(fastspec_spectrum_fl)
    assert isinstance(spec.data, HDF5RawSpectrum)
    check_obj(spec.data)


def test_read_acq(datadir):
    spec = Spectrum(datadir / "sample.acq")
    assert spec.data["spectra"]["Q"].shape == (32768, 1)
    assert spec.data["spectra"]["Q"].shape == spec.data["spectra"]["p0"].shape
    assert spec.data["spectra"]["Q"].shape == spec.data["spectra"]["p1"].shape
    assert spec.data["spectra"]["Q"].shape == spec.data["spectra"]["p2"].shape


def test_read_bad_format(datadir):
    spec = Spectrum(datadir / "bad.file")
    with pytest.raises(IOError):
        spec.data
