from os import path

import py
from edges_io import utils


def test_active_files(tmpdir: py.path.local):
    direc = tmpdir.mkdir("test_active_files")

    file1 = direc.join("this.txt")
    file2 = direc.join("that.txt")
    file3 = direc.join("ignored.old")
    file4 = direc.join("Notes.txt")

    file1.write("hey")
    file2.write("hey")
    file3.write("hey")
    file4.write("hey")

    fls = utils.get_active_files(direc.strpath)
    assert len(fls) == 2


def test_get_parent(tmpdir: py.path.local):
    direc = tmpdir.mkdir("test_get_parent").join("child").join("double_child")

    parent = utils.get_parent_dir(direc.strpath)
    assert path.basename(parent) == "child"
    root = utils.get_parent_dir(direc.strpath, 2)
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
