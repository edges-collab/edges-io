"""IO routines for VNA readings (S11, S12, S22 etc)."""

from __future__ import annotations

import warnings
from pathlib import Path

import attrs
import numpy as np
from astropy import units as un
from astropy.table import QTable

from . import types as tp


def _get_s1p_kind(path: Path) -> tuple[np.ndarray, str]:
    # identifying the format

    with path.open("r") as d:
        comment_rows = 0
        uses_commas = False
        flag = None
        lines = d.readlines()
        if lines[0].startswith("BEGIN") and lines[1].strip() in ["DB", "MA", "RI"]:
            # This is a format that has a BEGIN line, then a FLAG line, then data then END
            flag = lines[1].strip()
            comment_rows = 2
            footer_lines = 1
        else:
            for line in lines:
                # checking settings line
                if line.startswith("#"):
                    if "DB" in line or "dB" in line:
                        flag = "DB"
                    if "MA" in line:
                        flag = "MA"
                    if "RI" in line:
                        flag = "RI"

                    comment_rows += 1
                elif line.startswith(("!", "BEGIN")):
                    comment_rows += 1
                elif flag is not None:
                    if "," in line:
                        uses_commas = True
                    break
                else:
                    warnings.warn(
                        f"Non standard line in S11 file {path}: '{line}'\n"
                        "...Treating as a comment line.",
                        stacklevel=1,
                    )
                    comment_rows += 1

            # Also check the the last lines for stupid entries like "END"
            footer_lines = 0
            for line in lines[::-1]:
                if line.startswith("#") or ("END" in line) or not line:
                    footer_lines += 1
                else:
                    break

    if flag is None:
        raise OSError(f"The file {path} has incorrect format.")

    #  loading data
    d = np.genfromtxt(
        path,
        skip_header=comment_rows,
        skip_footer=footer_lines,
        delimiter="," if uses_commas else None,
    )

    return d, flag


def read_s1p(path: str | Path) -> QTable:
    """Read a file in either s1p or s2p format, recorded by a VNA.

    Parameters
    ----------
    path
        The path to the file.
    """
    d, flag = _get_s1p_kind(path)

    f = d[:, 0]
    table = QTable(
        {
            "frequency": f * un.Hz,
        }
    )

    if flag == "DB":
        table["s11"] = 10 ** (d[:, 1] / 20) * (
            np.cos((np.pi / 180) * d[:, 2]) + 1j * np.sin((np.pi / 180) * d[:, 2])
        )
    elif flag == "MA":
        table["s11"] = d[:, 1] * (
            np.cos((np.pi / 180) * d[:, 2]) + 1j * np.sin((np.pi / 180) * d[:, 2])
        )
    elif flag == "RI":
        table["s11"] = d[:, 1] + 1j * d[:, 2]

        if d.shape[1] > 3:
            table["s12"] = d[:, 3] + 1j * d[:, 4]
            table["s21"] = d[:, 5] + 1j * d[:, 6]
            table["s22"] = d[:, 7] + 1j * d[:, 8]

    else:
        raise ValueError("file had no flags set!")

    return table


@attrs.define
class SParams:
    freq: tp.FreqType = attrs.field()
    s11: np.ndarray | None = attrs.field()
    s12: np.ndarray | None = attrs.field()
    s21: np.ndarray | None = attrs.field()
    s22: np.ndarray | None = attrs.field()

    @s11.validator
    @s12.validator
    @s21.validator
    @s22.validator
    def _check_dims(self, attribute, value):
        if value is not None and value.shape != self.freq.shape:
            raise ValueError(
                f"Shape of {attribute.name} does not match shape of frequency"
            )

    @freq.validator
    def _check_freq(self, attribute, value):
        if not value.unit.is_equivalent(un.Hz):
            raise ValueError("Frequency must be in units of Hz")

        if value.ndim != 1:
            raise ValueError("Frequency must be 1D")

    @classmethod
    def from_table(cls, table: QTable):
        """Create an SParams object from a table."""
        # We slice each entry here so that we copy values, so we don't use a weakref
        return cls(
            freq=table["frequency"][:],
            s11=table["s11"][:] if "s11" in table.columns else None,
            s12=table["s12"][:] if "s12" in table.columns else None,
            s21=table["s21"][:] if "s21" in table.columns else None,
            s22=table["s22"][:] if "s22" in table.columns else None,
        )

    def to_table(self) -> QTable:
        """Convert to a table."""
        return QTable(self.to_dict())

    def to_dict(self) -> dict[str, np.ndarray]:
        """Convert to a dictionary."""
        return {k: v for k, v in attrs.asdict(self).items() if v is not None}

    @classmethod
    def from_s1p_file(cls, path: str | Path):
        """Read an S1P file."""
        table = read_s1p(path)
        return cls.from_table(table)
