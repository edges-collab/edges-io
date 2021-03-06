import pytest

import numpy as np
from pathlib import Path

from edges_io.h5 import HDF5RawSpectrum
from edges_io.io import Spectrum


def test_hdf5rawspectrum(fastspec_spectrum_fl):
    obj = HDF5RawSpectrum(fastspec_spectrum_fl)

    assert "spectra" in obj
    assert "freq_ancillary" in obj
    assert "time_ancillary" in obj

    # Look at an item
    b = obj["spectra"]["Q"][...]
    assert b.ndim == 2
    assert "Q" not in obj["spectra"].__memcache__

    # Load an item
    Q = obj["spectra"].load("Q")
    assert np.all(Q == b)
    assert "Q" in obj["spectra"].__memcache__
    assert "p0" not in obj["spectra"].__memcache__

    # Now load everything
    spectra = obj.load("spectra")
    Q2 = spectra["Q"]
    assert np.all(Q2 == Q)
    assert "p0" in obj["spectra"].__memcache__
    assert "time_ancillary" not in obj.__memcache__


def test_io_read(fastspec_spectrum_fl):
    spec = Spectrum(fastspec_spectrum_fl)

    assert isinstance(spec.data, HDF5RawSpectrum)

    obj = spec.data

    assert "spectra" in obj
    assert "freq_ancillary" in obj
    assert "time_ancillary" in obj

    # Look at an item
    b = obj["spectra"]["Q"][...]
    assert b.ndim == 2
    assert "Q" not in obj["spectra"].__memcache__

    # Load an item
    Q = obj["spectra"].load("Q")
    assert np.all(Q == b)
    assert "Q" in obj["spectra"].__memcache__
    assert "p0" not in obj["spectra"].__memcache__

    # Now load everything
    spectra = obj.load("spectra")
    Q2 = spectra["Q"]
    assert np.all(Q2 == Q)
    assert "p0" in obj["spectra"].__memcache__
    assert "time_ancillary" not in obj.__memcache__


def test_read_acq(datadir):
    spec = Spectrum(datadir / "sample.acq")
    assert spec.data["spectra"]["Q"].shape == (32768, 1)
    assert spec.data["spectra"]["Q"].shape == spec.data["spectra"]["p0"].shape
    assert spec.data["spectra"]["Q"].shape == spec.data["spectra"]["p1"].shape
    assert spec.data["spectra"]["Q"].shape == spec.data["spectra"]["p2"].shape


def test_read_bad_format(datadir):
    spec = Spectrum(datadir / "bad.file")
    with pytest.raises(ValueError):
        spec.data
