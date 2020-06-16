# -*- coding: utf-8 -*-
from pathlib import Path

import pytest

import h5py
import numpy as np
from edges_io.h5 import HDF5RawSpectrum


@pytest.fixture(scope="session")
def tmpdir(tmp_path_factory):
    return tmp_path_factory.mktemp("edges-io-tests")


@pytest.fixture(scope="session", autouse=True)
def fastspec_spectrum_fl(tmpdir):
    """An auto-generate empty Fastspec h5 format file"""

    ntimes = 2
    nfreqs = 32768
    flname = tmpdir / "fastspec_example_file.h5"

    attrs = {}
    attrs["fastspec_version"] = "0.0.0"
    attrs["start"] = 0
    attrs["stop"] = 0
    attrs["site"] = "simulated"
    attrs["instrument"] = "mid"
    attrs["switch_io_port"] = 57360
    attrs["switch_delay"] = 0.5
    attrs["input_channel"] = 1
    attrs["voltage_range"] = 0
    attrs["samples_per_accumulation"] = 4294967269
    attrs["acquisition_rate"] = 400
    attrs["num_channels"] = nfreqs
    attrs["num_taps"] = 5
    attrs["window_function_id"] = 3
    attrs["num_fft_threads"] = 10
    attrs["num_fft_buffers"] = 100
    attrs["stop_cycles"] = 0
    attrs["stop_seconds"] = 0.0
    attrs["stop_time"] = ""

    spec = {}
    spec["p0"] = np.zeros((ntimes, nfreqs))
    spec["p1"] = np.zeros((ntimes, nfreqs))
    spec["p2"] = np.zeros((ntimes, nfreqs))
    spec["Q"] = np.zeros((ntimes, nfreqs))

    fq_anc = {}
    fq_anc["frequencies"] = np.linspace(40, 200, nfreqs)

    time_anc = {}
    time_anc["times"] = np.array(
        ["2020:001:01:01:01", "2020:001:01:02:01"], dtype="|S17"
    )
    time_anc["adcmax"] = np.zeros((ntimes, 3))
    time_anc["adcmin"] = np.zeros((ntimes, 3))
    time_anc["data_drops"] = np.zeros((ntimes, 3), dtype="int")

    spectrum = HDF5RawSpectrum.from_data(
        {
            "meta": attrs,
            "spectra": spec,
            "time_ancillary": time_anc,
            "freq_ancillary": fq_anc,
        }
    )

    spectrum.write(flname)
    return flname


@pytest.fixture(scope="session")
def datadir() -> Path:
    return Path(__file__).parent / "test_data"
