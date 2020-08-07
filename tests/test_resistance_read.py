import shutil
from pathlib import Path

from edges_io.io import Resistance


def test_resistance_read_header(datadir: Path, tmpdir: Path):

    header = Resistance.read_old_style_csv_header(datadir / "resistance_header.csv")

    assert header["Start Time"] == "9/17/2017 10:25:01 AM"

    shutil.copyfile(
        datadir / "resistance_header.csv", tmpdir / "Ambient_1_2019_150_lab.csv"
    )

    path, _ = Resistance.check_self(tmpdir / "Ambient_1_2019_150_lab.csv", fix=True)

    assert path.name == "Ambient_01_2017_260_10_25_01_lab.csv"
