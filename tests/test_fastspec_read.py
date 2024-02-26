from __future__ import annotations

import pytest
from edges_io.io import Spectrum
from pygsdata import GSData
from read_acq.gsdata import read_acq_to_gsdata


def test_hdf5rawspectrum(fastspec_spectrum_fl, edgesloc):
    obj = read_acq_to_gsdata(fastspec_spectrum_fl, telescope_location=edgesloc)
    assert isinstance(obj, GSData)


def test_io_read(fastspec_spectrum_fl):
    spec = Spectrum(fastspec_spectrum_fl)
    assert isinstance(spec.get_data(), GSData)


def test_read_bad_format(datadir):
    spec = Spectrum(datadir / "bad.file")
    with pytest.raises(OSError, match="does not exist"):
        spec.get_data()
