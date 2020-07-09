import py
from os import path
from pathlib import Path

from edges_io import utils


def test_active_files(tmpdir: Path):
    direc = tmpdir / "test_active_files"
    direc.mkdir()

    file1 = direc / "this.txt"
    file2 = direc / "that.txt"
    file3 = direc / "ignored.old"
    file4 = direc / "Notes.txt"

    file1.touch()
    file2.touch()
    file3.touch()
    file4.touch()

    fls = utils.get_active_files(direc)
    assert len(fls) == 2


def test_get_parent(tmpdir: Path):
    direc = tmpdir / "test_get_parent/child/double_child"
    direc.mkdir(parents=True)

    parent = utils.get_parent_dir(str(direc))
    assert path.basename(parent) == "child"
    root = utils.get_parent_dir(str(direc), 2)
    assert path.basename(root) == "test_get_parent"


def test_ymd_to_jd():
    jd = utils.ymd_to_jd(2019, 1, 1)
    assert jd == 1

    jd = utils.ymd_to_jd(2019, 1, 30)
    assert jd == 30

    jd = utils.ymd_to_jd(2019, 3, 1)
    assert jd == 60

    # Ensure leap years go correctly
    jd = utils.ymd_to_jd(2020, 3, 1)
    assert jd == 61
