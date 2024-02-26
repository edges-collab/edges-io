import copy
from pathlib import Path

import numpy as np
import pytest
from astropy import units as apu
from astropy.coordinates import EarthLocation
from astropy.time import Time
from edges_io import TEST_DATA_PATH
from pygsdata import GSData
from read_acq.gsdata import write_gsdata_to_acq


@pytest.fixture(scope="session")
def tmpdir(tmp_path_factory):
    return tmp_path_factory.mktemp("edges-io-tests")


@pytest.fixture(scope="session")
def edgesloc():
    return EarthLocation(lat=-26.714778 * apu.deg, lon=116.605528 * apu.deg)


@pytest.fixture(scope="session")
def small_gsdata_obj(edgesloc: EarthLocation):
    ntimes = 2
    nfreqs = 32768
    npols = 1
    nloads = 3
    return GSData(
        data=np.zeros((nloads, npols, ntimes, nfreqs)),
        freq_array=np.linspace(0, 200, nfreqs) * apu.MHz,
        time_array=Time(
            [
                ["2020:001:01:01:01", "2020:001:01:01:01", "2020:001:01:01:01"],
                ["2020:001:01:02:01", "2020:001:01:02:01", "2020:001:01:02:01"],
            ]
        ),
        telescope_location=edgesloc,
        data_unit="power",
        auxiliary_measurements={
            "adcmax": np.zeros((ntimes, nloads)),
            "adcmin": np.zeros((ntimes, nloads)),
            "data_drops": np.zeros((ntimes, nloads), dtype="int"),
        },
    )


@pytest.fixture(scope="session", autouse=True)
def fastspec_spectrum_fl(tmpdir, small_gsdata_obj: GSData):
    """An auto-generated empty Fastspec h5 format file."""
    flname = tmpdir / "fastspec_example_file.acq"
    write_gsdata_to_acq(small_gsdata_obj, flname)
    return flname


@pytest.fixture(scope="session")
def datadir() -> Path:
    return TEST_DATA_PATH
