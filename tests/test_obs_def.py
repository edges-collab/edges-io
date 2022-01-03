"""Tests for defining observations inside the observation itself."""
import os
from pathlib import Path

from edges_io.io import CalibrationObservation


def test_trivial(datadir: Path):
    """Test reading in an already fully defined observation."""
    calobs = CalibrationObservation(datadir / "Receiver01_25C_2019_11_26_040_to_200MHz")
    assert calobs.path.name == "Receiver01_25C_2019_11_26_040_to_200MHz"


def test_include(datadir: Path):
    """Test reading in an observation which includes most of another one."""
    calobs = CalibrationObservation.from_def(
        datadir / "Receiver01_25C_2020_11_26_040_to_200MHz"
    )

    assert "Receiver01_25C_2020_11_26_040_to_200MHz" in os.listdir(calobs._tmpdir)

    assert "Receiver01_25C_2020_11_26_040_to_200MHz" in str(
        (calobs.path / "Resistance/AntSim4_01_2019_338_16_55_22_lab.csv").resolve()
    )
    assert "Receiver01_25C_2020_11_26_040_to_200MHz" in str(
        (calobs.path / "S11/AntSim401/External01.s1p").resolve()
    )

    # It's easy to not get the upper run num because they have the same type.
    # Let's make sure we get it.
    assert (calobs.path / "S11/AntSim401/External02.s1p").exists()

    # Ensure included files are from other observation.
    assert "Receiver01_25C_2019_11_26_040_to_200MHz" in str(
        (calobs.path / "S11/AntSim301/External01.s1p").resolve()
    )


def test_prefer(datadir: Path):
    """Test reading in an observation which _prefers_ another one."""
    calobs = CalibrationObservation.from_def(
        datadir / "Receiver01_25C_2021_11_26_040_to_200MHz", include_previous=False
    )

    assert "Receiver01_25C_2021_11_26_040_to_200MHz" in os.listdir(calobs._tmpdir)

    assert "Receiver01_25C_2020_11_26_040_to_200MHz" in str(
        (calobs.path / "Resistance/AntSim4_01_2019_338_16_55_22_lab.csv").resolve()
    )
    assert "Receiver01_25C_2020_11_26_040_to_200MHz" in str(
        (calobs.path / "S11/AntSim401/External01.s1p").resolve()
    )

    # It's easy to not get the upper run num because they have the same type.
    # Let's make sure we get it.
    assert (calobs.path / "S11/AntSim401/External02.s1p").exists()

    # Ensure included files are from other observation.
    assert "Receiver01_25C_2021_11_26_040_to_200MHz" in str(
        (calobs.path / "S11/AntSim301/External01.s1p").resolve()
    )


def test_default_include(datadir: Path):
    calobs = CalibrationObservation.from_def(
        datadir / "Receiver01_25C_2019_12_26_040_to_200MHz"
    )

    assert "Receiver01_25C_2019_12_26_040_to_200MHz" in os.listdir(calobs._tmpdir)

    assert "Receiver01_25C_2019_12_26_040_to_200MHz" in str(
        (calobs.path / "Resistance/AntSim4_01_2019_338_16_55_22_lab.csv").resolve()
    )
    assert "Receiver01_25C_2019_12_26_040_to_200MHz" in str(
        (calobs.path / "S11/AntSim401/External01.s1p").resolve()
    )

    # It's easy to not get the upper run num because they have the same type.
    # Let's make sure we get it.
    assert (calobs.path / "S11/AntSim401/External02.s1p").exists()

    # Ensure included files are from other observation.
    assert "Receiver01_25C_2019_11_26_040_to_200MHz" in str(
        (calobs.path / "S11/AntSim301/External01.s1p").resolve()
    )


def test_external_observation(datadir: Path):
    calobs = CalibrationObservation.from_observation_yaml(datadir / "observation.yaml")

    assert (calobs.path / "S11/HotLoad01/External02.s1p").exists()

    assert "Receiver01_25C_2019_11_26_040_to_200MHz" in str(
        (calobs.path / "S11/HotLoad01/External02.s1p").resolve()
    )

    assert "Receiver01_25C_2020_11_26_040_to_200MHz" in str(
        (calobs.path / "S11/SwitchingState01/ExternalMatch02.s1p").resolve()
    )
