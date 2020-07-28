from pathlib import Path

from edges_io.io import Spectrum


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
