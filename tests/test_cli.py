import pytest

import shutil
from click.testing import CliRunner
from pathlib import Path
from typing import Tuple

from edges_io import cli
from edges_io.utils import console


@pytest.mark.parametrize(
    "folder",
    [
        "Receiver01_25C_2019_11_26_040_to_200MHz",
        "Receiver01_25C_2019_12_26_040_to_200MHz",
        "Receiver01_25C_2020_11_26_040_to_200MHz",
        "Receiver01_25C_2021_11_26_040_to_200MHz",
    ],
)
def test_check(datadir, caplog, folder):
    runner = CliRunner()
    result = runner.invoke(cli.check, [str(datadir / folder)])

    assert result.exit_code == 0
    assert "SUCCESS" in caplog.records[-1].levelname


@pytest.mark.parametrize("fix_strategy", ["y", "i", "o", "p"])
def test_check_fix(datadir, tmpdir, caplog, monkeypatch, fix_strategy):
    runner = CliRunner()
    folder = "Receiver01_25C_2022_11_26_040_to_200MHz"
    tmpdir = tmpdir / fix_strategy
    tmpdir.mkdir()
    shutil.copytree(datadir / folder, tmpdir / folder)

    # Patch the console.input() function to always return 'y' -- which will delete
    # stuff in the _ask_to_rm function.
    monkeypatch.setattr(
        console,
        "input",
        lambda msg: "badfile.here.old" if "Change " in msg else fix_strategy,
    )

    result = runner.invoke(cli.check, [str(tmpdir / folder), "--fix"])

    assert result.exit_code == 0
    assert "SUCCESS" in caplog.records[-1].levelname


def test_check_verbosity_noop(datadir, caplog):
    runner = CliRunner()

    result = runner.invoke(
        cli.check, [str(datadir / "Receiver01_25C_2019_11_26_040_to_200MHz")]
    )

    assert result.exit_code == 0

    txt = caplog.text
    n = len(txt)

    # This adds and subtracts verbosity
    result = runner.invoke(
        cli.check, [str(datadir / "Receiver01_25C_2019_11_26_040_to_200MHz"), "-vV"]
    )
    assert result.exit_code == 0

    assert caplog.text[n:] == txt


def test_check_verbosity_extra(datadir, caplog):
    runner = CliRunner()

    result = runner.invoke(
        cli.check, [str(datadir / "Receiver01_25C_2019_11_26_040_to_200MHz")]
    )

    assert result.exit_code == 0

    txt = caplog.text
    n = len(txt)

    # This subtracts verbosity
    result = runner.invoke(
        cli.check, [str(datadir / "Receiver01_25C_2019_11_26_040_to_200MHz"), "-VVV"]
    )
    assert result.exit_code == 0

    assert caplog.text[n:] != txt


def test_check_verbosity_overkill(datadir, caplog):
    runner = CliRunner()

    result = runner.invoke(
        cli.check, [str(datadir / "Receiver01_25C_2019_11_26_040_to_200MHz"), "-vvvv"]
    )

    assert result.exit_code == 0

    # This subtracts verbosity
    result = runner.invoke(
        cli.check,
        [str(datadir / "Receiver01_25C_2019_11_26_040_to_200MHz"), "-VVVVVVV"],
    )
    assert result.exit_code == 0


def unmove(temp: str, datadir: Path, tmpdir: Path) -> Tuple[str, Path]:
    folder = "Receiver01_25C_2019_11_26_040_to_200MHz"
    bad = tmpdir / "Receiver01_2019_11_26_040_to_200MHz"
    shutil.copytree(datadir / folder, bad / temp)
    return folder, bad


def test_mv(datadir: Path, tmpdir: Path):
    folder, bad = unmove("25C", datadir, tmpdir)

    runner = CliRunner()
    result = runner.invoke(cli.mv, [str(bad / "25C"), "--clean"])

    assert result.exit_code == 0
    assert not bad.exists()
    assert (tmpdir / folder).exists()
    assert not (tmpdir / f"{folder}/25C").exists()
    assert (tmpdir / f"{folder}/Spectra").exists()


def test_mv_all(datadir: Path, tmpdir: Path):
    tmpdir = tmpdir / "mvall"
    folder, bad = unmove("25C", datadir, tmpdir)
    unmove("35C", datadir, tmpdir)
    unmove("15C", datadir, tmpdir)

    print(list((bad / "25C").glob("*")))

    runner = CliRunner()
    folders = [str(x) for x in bad.glob("*C")]
    result = runner.invoke(cli.mv_all, folders + ["--clean"])
    print(result.stdout)
    print(result.exception)

    assert result.exit_code == 0
    assert not bad.exists()
    assert (tmpdir / folder).exists()
    assert not (tmpdir / f"{folder}/25C").exists()
    assert not (tmpdir / f"{folder}/35C").exists()
    assert not (tmpdir / f"{folder}/15C").exists()

    assert (tmpdir / f"{folder}/Spectra").exists()
