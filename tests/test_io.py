import pytest

import logging
from bidict import bidict
from pathlib import Path

from edges_io import io, utils

LOAD_ALIASES = bidict(
    {
        "ambient": "Ambient",
        "hot_load": "HotLoad",
        "open": "LongCableOpen",
        "short": "LongCableShorted",
    }
)

LOGGING = logging.getLogger("edges-io")


@pytest.fixture(scope="module")
def test_dir(tmp_path_factory):
    return test_env(tmp_path_factory)


@pytest.fixture(scope="module")
def test_env(tmp_path_factory):
    # Create an ideal observation file using tmp_path_factory
    pthList = ["Spectra", "Resistance", "S11"]
    s11List = [
        "Ambient01",
        "AntSim301",
        "HotLoad01",
        "LongCableOpen01",
        "LongCableShorted01",
        "ReceiverReading01",
        "ReceiverReading02",
        "SwitchingState01",
        "SwitchingState02",
    ]
    root_dir = tmp_path_factory.mktemp("Test_Obs")
    obs_dir = root_dir / "Receiver01_25C_2020_01_01_010_to_200MHz"
    obs_dir.mkdir()
    note = obs_dir / "Notes.txt"
    note.touch()
    dlist = []
    slist = []
    for i, p in enumerate(pthList):
        dlist.append(obs_dir / p)
        dlist[i].mkdir()
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
                name1 = filename + "_01_2020_001_01_01_01_lab.csv"
                file1 = dlist[i] / name1
                file1.touch()
        elif p == "S11":
            print("Making S11 files")
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
                    name1 = filename + "01.s1p"
                    name2 = filename + "02.s1p"
                    file1 = slist[k] / name1
                    file1.write_text(
                        "# Hz S RI R 50\n"
                        "40000000        0.239144887761343       0.934085904901478\n"
                        "40000000        0.239144887761343       0.934085904901478"
                    )
                    file2 = slist[k] / name2
                    file2.write_text(
                        "# Hz S RI R 50\n"
                        "40000000        0.239144887761343       0.934085904901478\n"
                        "40000000        0.239144887761343       0.934085904901478"
                    )

        elif p == "Spectra":
            print("Making Spectra files")
            fileList = [
                "Ambient",
                "AntSim3",
                "HotLoad",
                "LongCableOpen",
                "LongCableShorted",
            ]
            for filename in fileList:
                name1 = filename + "_01_2020_001_01_01_01_lab.acq"
                file1 = dlist[i] / name1
                file1.touch()
    return obs_dir


# function to make observation object
def new_test_obs(testdir):
    return io.CalibrationObservation(
        testdir, include_previous=False, compile_from_def=False
    )


# directory testing
def test_make_good_obs(test_env, caplog):
    # test that correct layouts pass (make an obs)
    new_test_obs(test_env)


def test_bad_dirname_obs(test_env, caplog):
    # test that incorrect directories fail
    test_dir = test_env
    base = test_dir.parent
    wrong_dir = base / "Receiver_2020_01_01_010_to_200MHz"
    test_dir.rename(wrong_dir)
    with pytest.raises(utils.FileStructureError):
        new_test_obs(wrong_dir)
    print(caplog.text)
    assert (
        "The filename Receiver_2020_01_01_010_to_200MHz does not have the correct format"
        in caplog.text
    )

    # receiver number
    test_dir = wrong_dir
    wrong_dir = base / "Receiver00_25C_2020_01_01_010_to_200MHz"
    test_dir.rename(wrong_dir)
    print("WRONGDIR: ", wrong_dir)
    new_test_obs(wrong_dir)
    assert "Unknown receiver number" in caplog.text

    # year
    test_dir = wrong_dir
    wrong_dir = base / "Receiver01_25C_2009_01_01_010_to_200MHz"
    test_dir.rename(wrong_dir)
    new_test_obs(wrong_dir)
    assert "Unknown year" in caplog.text

    test_dir = wrong_dir
    wrong_dir = base / "Receiver01_25C_2045_01_01_010_to_200MHz"
    test_dir.rename(wrong_dir)
    new_test_obs(wrong_dir)
    assert "Unknown year" in caplog.text

    # month
    test_dir = wrong_dir
    wrong_dir = base / "Receiver01_25C_2020_13_01_010_to_200MHz"
    test_dir.rename(wrong_dir)
    new_test_obs(wrong_dir)
    assert "Unknown month" in caplog.text

    # day
    test_dir = wrong_dir
    wrong_dir = base / "Receiver01_25C_2020_01_32_010_to_200MHz"
    test_dir.rename(wrong_dir)
    new_test_obs(wrong_dir)
    assert "Unknown day" in caplog.text

    # freqlow
    test_dir = wrong_dir
    wrong_dir = base / "Receiver01_25C_2020_01_01_000_to_200MHz"
    test_dir.rename(wrong_dir)
    new_test_obs(wrong_dir)
    assert "Low frequency is weird" in caplog.text

    # freqhigh
    test_dir = wrong_dir
    wrong_dir = base / "Receiver01_25C_2020_01_01_010_to_900MHz"
    test_dir.rename(wrong_dir)
    new_test_obs(wrong_dir)
    assert "High frequency is weird" in caplog.text

    # freqrange
    test_dir = wrong_dir
    wrong_dir = base / "Receiver01_25C_2020_01_01_200_to_010MHz"
    test_dir.rename(wrong_dir)
    new_test_obs(wrong_dir)
    assert "Low frequency > High Frequency" in caplog.text


def test_spectra_run_num(datadir: Path):
    spec = io.Spectra(
        datadir / "Receiver01_25C_2019_11_26_040_to_200MHz/Spectra", run_num=1
    )
    assert isinstance(spec.run_num, dict)
    assert all(int(v) == 1 for v in spec.run_num.values())


def test_resistance_run_num(datadir: Path):
    spec = io.Resistances(
        datadir / "Receiver01_25C_2019_11_26_040_to_200MHz/Resistance", run_num=1
    )
    assert isinstance(spec.run_num, dict)
    assert all(int(v) == 1 for v in spec.run_num.values())


def test_list_of_files(datadir: Path):
    obs = datadir / "Receiver01_25C_2019_11_26_040_to_200MHz"
    calobs = io.CalibrationObservation(obs)

    lof = [fl.relative_to(obs.parent) for fl in calobs.list_of_files]

    for fl in lof:
        print(fl)
    assert (
        obs.relative_to(obs.parent)
        / "Resistance"
        / "Ambient_01_2019_329_16_02_35_lab.csv"
        in lof
    )
    assert obs.relative_to(obs.parent) / "S11" / "Ambient01" / "External02.s1p" in lof
    assert (
        obs.relative_to(obs.parent) / "S11" / "Ambient02" / "External02.s1p" not in lof
    )
    assert (
        obs.relative_to(obs.parent) / "S11" / "Ambient01" / "External01.s1p" not in lof
    )


### read testing
# do identical acq and h5 files read in identically?
# how does the read handle missing fields?


### spectra testing
# test to see if frequency range matches spectra shape

# test to see if flow, fhigh match file

### resistance testing
