import pytest

import h5py
import numpy as np
from pathlib import Path

from edges_io.h5 import HDF5RawSpectrum


@pytest.fixture(scope="session")
def tmpdir(tmp_path_factory):
    return tmp_path_factory.mktemp("edges-io-tests")


@pytest.fixture(scope="session")
def fastspec_data():
    ntimes = 2
    nfreqs = 32768

    attrs = {
        "fastspec_version": "0.0.0",
        "start": 0,
        "stop": 0,
        "site": "simulated",
        "instrument": "mid",
        "switch_delay": 0.5,
        "input_channel": 1,
        "voltage_range": 0,
        "samples_per_accumulation": 4294967269,
        "acquisition_rate": 400,
        "num_channels": nfreqs,
        "num_taps": 5,
        "window_function_id": 3,
        "num_fft_threads": 10,
        "num_fft_buffers": 100,
        "stop_cycles": 0,
        "stop_seconds": 0.0,
        "stop_time": "",
    }

    spec = {
        "p0": np.zeros((ntimes, nfreqs)),
        "p1": np.zeros((ntimes, nfreqs)),
        "p2": np.zeros((ntimes, nfreqs)),
        "Q": np.zeros((ntimes, nfreqs)),
    }

    fq_anc = {"frequencies": np.linspace(40, 200, nfreqs)}
    time_anc = {
        "times": np.array(["2020:001:01:01:01", "2020:001:01:02:01"], dtype="|S17"),
        "adcmax": np.zeros((ntimes, 3)),
        "adcmin": np.zeros((ntimes, 3)),
        "data_drops": np.zeros((ntimes, 3), dtype="int"),
    }

    return {
        "meta": attrs,
        "spectra": spec,
        "time_ancillary": time_anc,
        "freq_ancillary": fq_anc,
    }


@pytest.fixture(scope="session", autouse=True)
def fastspec_spectrum_fl(tmpdir, fastspec_data):
    """An auto-generate empty Fastspec h5 format file"""
    spectrum = HDF5RawSpectrum.from_data(fastspec_data)
    flname = tmpdir / "fastspec_example_file.h5"

    spectrum.write(flname)
    return flname


@pytest.fixture(scope="session")
def datadir() -> Path:
    return Path(__file__).parent / "test_data"
