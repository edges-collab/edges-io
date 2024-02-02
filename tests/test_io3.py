import pytest

from astropy import units as un
from pathlib import Path

from edges_io import io3
from edges_io.vna import read_s1p


@pytest.fixture(scope="module")
def smallcal(datadir: Path) -> io3.CalibrationObservation:
    return io3.CalibrationObservation.from_date(
        root_dir=datadir / "edges3-mock-root",
        year=2023,
        day=70,
    )


def test_all_files_present(smallcal: io3.CalibrationObservation):
    loads = ["hot", "amb", "short", "open"]
    assert len(smallcal.s11_files) == len(loads) + 1  # +1 for the LNA
    assert all(load in smallcal.s11_files for load in loads)

    assert all(
        p.exists()
        for loadfiles in smallcal.s11_files.values()
        for p in loadfiles.values()
    )

    assert all(load in smallcal.acq_files for load in loads)
    assert all(f.exists() for f in smallcal.acq_files.values())
    assert smallcal.temperature_file.exists()


def test_temperature_read(smallcal: io3.CalibrationObservation):
    temp_table = smallcal.get_temperature_table()

    assert "time" in temp_table.columns
    assert "hot_load_temperature" in temp_table.columns
    assert "amb_load_temperature" in temp_table.columns
    assert "front_end_temperature" in temp_table.columns
    assert "inner_box_temperature" in temp_table.columns
    assert "thermal_control" in temp_table.columns
    assert "battery_voltage" in temp_table.columns
    assert "pr59_current" in temp_table.columns

    for col in temp_table.columns:
        if col.endswith("temperature"):
            assert temp_table[col].unit == un.K


def test_get_mean_temperature(smallcal: io3.CalibrationObservation):
    temp_table = smallcal.get_temperature_table()
    mean_temp = io3.get_mean_temperature(temp_table, load="amb")
    assert mean_temp.unit == un.K
    assert mean_temp.value == pytest.approx(300, abs=15)

    mean_temp0 = io3.get_mean_temperature(
        temp_table, load="amb", start_time=temp_table["time"].min()
    )
    assert mean_temp0 == mean_temp

    mean_temp0 = io3.get_mean_temperature(
        temp_table, load="open", end_time=temp_table["time"].max()
    )
    assert mean_temp0 == mean_temp

    mean_temp0 = io3.get_mean_temperature(temp_table, load="hot")
    assert mean_temp0 > mean_temp

    with pytest.raises(ValueError):
        io3.get_mean_temperature(temp_table, load="fake")


def test_read_s11s(smallcal: io3.CalibrationObservation):
    files = sum((list(d.values()) for d in smallcal.s11_files.values()), start=[])
    for fl in files:
        table = read_s1p(fl)
        assert table["frequency"].unit == un.Hz
