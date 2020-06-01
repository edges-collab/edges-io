# -*- coding: utf-8 -*-
"""
    Dummy conftest.py for edges_io.

    If you don't know what this is for, just leave it empty.
    Read more about conftest.py under:
    https://pytest.org/latest/plugins.html
"""

import pytest

import h5py
import numpy as np


@pytest.fixture(scope="session")
def tmpdir(tmp_path_factory):
    return tmp_path_factory.mktemp("edges-io-tests")


@pytest.fixture(scope="session", autouse=True)
def fastspec_spectrum_fl(tmpdir):
    """An auto-generate empty Fastspec h5 format file"""

    ntimes = 2
    nfreqs = 32768

    flname = tmpdir / "fastspec_example_file.h5"
    with h5py.File(flname, "w") as fl:
        fl.attrs["fastspec_version"] = "0.0.0"
        fl.attrs["start"] = 0
        fl.attrs["stop"] = 0
        fl.attrs["site"] = "simulated"
        fl.attrs["instrument"] = "mid"
        fl.attrs["switch_io_port"] = 57360
        fl.attrs["switch_delay"] = 0.5
        fl.attrs["input_channel"] = 1
        fl.attrs["voltage_range"] = 0
        fl.attrs["samples_per_accumulation"] = 4294967269
        fl.attrs["acquisition_rate"] = 400
        fl.attrs["num_channels"] = nfreqs
        fl.attrs["num_taps"] = 5
        fl.attrs["window_function_id"] = 3
        fl.attrs["num_fft_threads"] = 10
        fl.attrs["num_fft_buffers"] = 100
        fl.attrs["stop_cycles"] = 0
        fl.attrs["stop_seconds"] = 0.0
        fl.attrs["stop_time"] = ""

        spec = fl.create_group("spectra")
        spec["p0"] = np.zeros((ntimes, nfreqs))
        spec["p1"] = np.zeros((ntimes, nfreqs))
        spec["p2"] = np.zeros((ntimes, nfreqs))
        spec["Q"] = np.zeros((ntimes, nfreqs))

        fq_anc = fl.create_group("freq_ancillary")
        fq_anc["frequencies"] = np.linspace(40, 200, nfreqs)

        time_anc = fl.create_group("time_ancillary")
        time_anc["times"] = np.array(
            ["2020:001:01:01:01", "2020:001:01:02:01"], dtype="|S17"
        )
        time_anc["adcmax"] = np.zeros((ntimes, 3))
        time_anc["adcmin"] = np.zeros((ntimes, 3))
        time_anc["data_drops"] = np.zeros((ntimes, 3), dtype="int")
    return flname
