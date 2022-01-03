import pytest

import numpy as np
import shutil
from pathlib import Path

from edges_io.io import Resistance


def test_resistance_read_old_header(datadir: Path, tmpdir: Path):

    header, nlines = Resistance.read_old_style_csv_header(
        datadir / "old_resistance_file.csv"
    )

    assert header["Start Time"] == "9/14/2017 2:16:45 PM"

    shutil.copyfile(
        datadir / "old_resistance_file.csv", tmpdir / "Ambient_1_2019_150_lab.csv"
    )

    path, _ = Resistance.check_self(tmpdir / "Ambient_1_2019_150_lab.csv", fix=True)

    assert path.name == "Ambient_01_2017_257_14_16_45_lab.csv"


def test_resistance_read_new(datadir: Path):
    fl = (
        datadir
        / "Receiver01_25C_2019_11_26_040_to_200MHz/Resistance/Ambient_01_2019_329_16_02_35_lab.csv"
    )

    r = Resistance(fl)
    r.read()
    assert len(r.resistance) == 9
    assert len(r.resistance.dtype.names) == 12
    assert len(r.ancillary) == 0


def test_resistance_read_old(datadir: Path):
    fl = datadir / "old_resistance_file.csv"

    r = Resistance(fl, check=False)
    r.read()
    assert len(r.resistance) == 11
    assert len(r.resistance.dtype.names) == 11
    assert len(r.ancillary) == 0
    assert not np.any(np.isnan(r.resistance["load_resistance"]))
