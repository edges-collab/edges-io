import configparser
import logging
import os
import shutil

import pytest

import py
from bidict import bidict
from edges_io import io

LOAD_ALIASES = bidict(
    {
        "ambient": "Ambient",
        "hot_load": "HotLoad",
        "open": "LongCableOpen",
        "short": "LongCableShorted",
    }
)

LOGGING = logging.getLogger("edges-io")


@pytest.fixture("module")
def test_dir(tmp_path_factory):
    testDir = test_env(tmp_path_factory)
    return testDir


@pytest.fixture("module")
def test_env(tmp_path_factory):
    # Create an ideal observation file using tmp_path_factory
    pthList = ["Spectra", "Resistance", "S11"]
    s11List = [
        "Ambient",
        "AntSim3",
        "HotLoad",
        "LongCableOpen",
        "LongCableShorted",
        "ReceiverReading01",
        "ReceiverReading02",
        "SwitchingState01",
        "SwitchingState02",
    ]
    root_dir = tmp_path_factory.mktemp("Test_Obs")
    obs_dir = root_dir / "Receiver01_2020_01_01_010_to_200MHz"
    obs_dir.mkdir()
    temp_dir = obs_dir / "25C"
    temp_dir.mkdir()
    note = temp_dir / "Notes.txt"
    note.touch()
    dlist = []
    slist = []
    for i, p in enumerate(pthList):
        dlist.append(temp_dir / p)
        dlist[i].mkdir()
        if p == "Spectra":
            print("Making Spectra files")
            fileList = [
                "Ambient",
                "AntSim3",
                "HotLoad",
                "LongCableOpen",
                "LongCableShorted",
            ]
            for filename in fileList:
                # print(filename)
                name1 = filename + "_01_2020_001_01_01_01_lab.acq"
                file1 = dlist[i] / name1
                # print(file1)
                file1.touch()
        if p == "Resistance":
            print("Making Resistance files")
            fileList = [
                "Ambient",
                "AntSim3",
                "HotLoad",
                "LongCableOpen",
                "LongCableShorted",
            ]
            for filename in fileList:
                # print(filename)
                name1 = filename + "_01_2020_001_01_01_01_lab.csv"
                file1 = dlist[i] / name1
                # print(file1)
                file1.touch()
        if p == "S11":
            print("Making S11 files")
            fileList = ["External", "Match", "Open", "Short"]
            for k, s in enumerate(s11List):
                slist.append(dlist[i] / s)
                slist[k].mkdir()
                if s[:-2] == "ReceiverReading":
                    fileList = ["ReceiverReading", "Match", "Open", "Short"]
                elif s[:-2] == "SwitchingState":
                    fileList = [
                        "ExternalOpen",
                        "ExternalMatch",
                        "ExternalShort",
                        "Match",
                        "Open",
                        "Short",
                    ]
                else:
                    fileList = ["External", "Match", "Open", "Short"]
                for filename in fileList:
                    # print(filename)
                    name1 = filename + "01.s1p"
                    name2 = filename + "02.s1p"
                    file1 = slist[k] / name1
                    file1.write_text(
                        "# Hz S RI R 50\n40000000        0.239144887761343       0.934085904901478\n40000000        0.239144887761343       0.934085904901478"
                    )
                    # print(file1.read_text())
                    file2 = slist[k] / name2
                    file2.write_text(
                        "# Hz S RI R 50\n40000000        0.239144887761343       0.934085904901478\n40000000        0.239144887761343       0.934085904901478"
                    )

    return obs_dir


# function to make observation object
def new_testObs(testdir):
    testObs = io.CalibrationObservation(testdir)
    return testObs


### directory testing
def test_make_good_Obs(test_env, caplog):
    # test that correct layouts pass (make an obs)
    testDir = test_env
    new_testObs(testDir)


def test_bad_dirname_Obs(test_env, caplog):
    # test that incorrect directories fail
    testDir = test_env
    base = testDir.parent
    wrongDir = base / "Receiver_2020_01_01_010_to_200MHz"
    testDir.rename(wrongDir)
    with pytest.raises(Exception):
        new_testObs(wrongDir)
    assert "directory name is in the wrong format" in caplog.text
    # receiver number
    testDir = wrongDir
    wrongDir = base / "Receiver00_2020_01_01_010_to_200MHz"
    testDir.rename(wrongDir)
    with pytest.raises(Exception):
        new_testObs(wrongDir)
    assert "Unknown receiver number" in caplog.text
    # year
    testDir = wrongDir
    wrongDir = base / "Receiver01_2009_01_01_010_to_200MHz"
    testDir.rename(wrongDir)
    with pytest.raises(Exception):
        new_testObs(wrongDir)
    assert "Unknown year" in caplog.text
    testDir = wrongDir
    wrongDir = base / "Receiver01_2045_01_01_010_to_200MHz"
    testDir.rename(wrongDir)
    with pytest.raises(Exception):
        new_testObs(wrongDir)
    assert "Unknown year" in caplog.text
    # month
    testDir = wrongDir
    wrongDir = base / "Receiver01_2020_13_01_010_to_200MHz"
    testDir.rename(wrongDir)
    with pytest.raises(Exception):
        new_testObs(wrongDir)
    assert "Unknown month" in caplog.text
    # day
    testDir = wrongDir
    wrongDir = base / "Receiver01_2020_01_32_010_to_200MHz"
    testDir.rename(wrongDir)
    with pytest.raises(Exception):
        new_testObs(wrongDir)
    assert "Unknown day" in caplog.text
    # freqlow

    testDir = wrongDir
    wrongDir = base / "Receiver01_2020_01_01_000_to_200MHz"
    testDir.rename(wrongDir)
    with pytest.raises(Exception):
        new_testObs(wrongDir)
    assert "Low frequency is weird" in caplog.text
    # freqhigh
    testDir = wrongDir
    wrongDir = base / "Receiver01_2020_01_01_010_to_900MHz"
    testDir.rename(wrongDir)
    with pytest.raises(Exception):
        new_testObs(wrongDir)
    assert "High frequency is weird" in caplog.text
    # freqrange
    testDir = wrongDir
    wrongDir = base / "Receiver01_2020_01_01_200_to_010MHz"
    testDir.rename(wrongDir)
    with pytest.raises(Exception):
        new_testObs(wrongDir)
    assert "Low frequency > High Frequency" in caplog.text


### read testing
# do identical acq and h5 files read in identically?
# how does the read handle missing fields?


### spectra testing
# test to see if frequency range matches spectra shape

# test to see if flow, fhigh match file

### resistance testing
