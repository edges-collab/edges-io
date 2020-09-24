import pytest

from pathlib import Path

from edges_io.io import CalibrationObservation, Spectrum
from edges_io.utils import FileStructureError


def test_spectrum_file_params(tmpdir: Path):
    """test that the Spectrum class can use default values for file keys"""
    fname = tmpdir / "AmbientLoad_25C_01_01_2017_01_01_01.acq"
    fname.touch()

    path, _ = Spectrum.check_self(fname, fix=True)

    assert path.name == "Ambient_01_2017_001_01_01_01_lab.acq"

    fname = tmpdir / "SimAnt3_1_2017_150_lab.acq"
    fname.touch()

    path, _ = Spectrum.check_self(fname, fix=True)

    assert path.name == "AntSim3_01_2017_150_00_00_00_lab.acq"


def test_spectrum_file_param_validation(tmpdir: Path, caplog):
    fname = tmpdir / "Ambient_00_2050_400_61_61_61_lab.cst"
    fname.touch()

    path, _ = Spectrum.check_self(fname, fix=True)

    assert caplog.text.count("ERROR") == 7


def test_obs_to_time(datadir):
    fl = datadir / "Receiver01_25C_2019_11_26_040_to_200MHz"
    dt = CalibrationObservation.path_to_datetime(fl)
    assert dt.year == 2019
    assert dt.month == 11

    with pytest.raises(FileStructureError):
        CalibrationObservation.path_to_datetime("herpy_derpy")
