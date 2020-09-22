import pytest

import os
import shutil
from click.testing import CliRunner
from pathlib import Path

from edges_io import cli


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


def test_mv(datadir: Path, tmpdir: Path):
    folder = "Receiver01_25C_2019_11_26_040_to_200MHz"
    bad = tmpdir / "Receiver01_2019_11_26_040_to_200MHz"

    shutil.copytree(datadir / folder, bad)
    (bad / "25C").mkdir()
    shutil.move(bad / "Spectra", bad / "25C/Spectra")
    shutil.move(bad / "Resistance", bad / "25C/Resistance")
    shutil.move(bad / "S11", bad / "25C/S11")
    shutil.move(bad / "Notes.txt", bad / "25C/Notes.txt")

    runner = CliRunner()
    result = runner.invoke(cli.mv, [str(bad / "25C"), "--clean"])

    assert result.exit_code == 0
    assert not bad.exists()
    assert (tmpdir / folder).exists()
    assert not (tmpdir / f"{folder}/25C").exists()
    assert (tmpdir / f"{folder}/Spectra").exists()
