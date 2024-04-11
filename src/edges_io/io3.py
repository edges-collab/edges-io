"""Methods for dealing with EDGES-3 files and structures."""
from __future__ import annotations

import re
import warnings
from datetime import datetime, timedelta, timezone
from functools import cache
from pathlib import Path
from typing import Any, Literal

import attrs
import numpy as np
from astropy import units as un
from astropy.table import QTable
from astropy.time import Time
from frozendict import frozendict
from hickleable import hickleable

from . import io


def _get_single_s1p_file(
    root: Path,
    year: int,
    day: int,
    label: str,
    hour: int | str = "first",
    allow_closest: int = 30,
) -> Path:
    if isinstance(hour, int):
        glob = f"{year:04d}_{day:03d}_{hour:02d}_{label}.s1p"
    else:
        glob = f"{year:04d}_{day:03d}_??_{label}.s1p"

    file_temp = sorted(root.glob(glob))

    if len(file_temp) == 0:
        if allow_closest:
            for dday in range(1, allow_closest):
                try:
                    fl = _get_single_s1p_file(
                        root,
                        year,
                        day + dday,
                        label=label,
                        hour=hour,
                        allow_closest=False,
                    )
                    break
                except FileNotFoundError:
                    pass

                try:
                    fl = _get_single_s1p_file(
                        root,
                        year,
                        day - dday,
                        label=label,
                        hour=hour,
                        allow_closest=False,
                    )
                    break
                except FileNotFoundError:
                    pass
            return fl
        raise FileNotFoundError(f"No s1p files found in {root} with glob {glob}")
    elif len(file_temp) > 1:
        if hour == "first":
            return file_temp[0]
        elif hour == "last":
            return file_temp[-1]
        else:
            raise OSError(f"More than one file found for {year}, {day}, {label}")

    return file_temp[0]


def get_s1p_files(
    load: Literal["amb", "ant", "hot", "open", "short", "lna"],
    year: int,
    day: int,
    root_dir: Path | str,
    hour: int | str = "first",
    allow_closest_s11_within: int = 5,
) -> dict[str, Path]:
    """Take the load and return a list of .s1p files for that load."""
    root_dir = Path(root_dir)

    files = {
        "input": _get_single_s1p_file(
            root_dir, year, day, load, hour, allow_closest_s11_within
        )
    }

    for name, label in {"open": "O", "short": "S", "match": "L"}.items():
        if load == "lna":
            label = f"{load}_{label}"
        files[name] = _get_single_s1p_file(
            root_dir, year, day, label, hour, allow_closest_s11_within
        )

    return files


def get_acq_file(
    load: Literal["amb", "hot", "short", "open"], root: Path, year: int, day: int
) -> Path:
    """Get the ACQ files for a particular load."""
    d = root / "mro" / load / str(year)
    glob = f"{year}_{day:03d}_??_??_??_{load}.acq"
    files = sorted(d.glob(glob))

    if not files:
        raise FileNotFoundError(f"No files found in {d} for {glob}")
    elif len(files) > 1:
        raise OSError(f"More than one file found in {d} for {glob}")

    return files[0]


def read_temperature_log_entry(lines: list[str]) -> dict[str, Any]:
    if len(lines) != 10:
        raise ValueError("Expected 10 lines for a temperature log entry")

    # first line is the date
    year, day, hour = (int(x) for x in lines[0].split("_"))

    # second line is the time
    _, _, _, timestamp, _, _year = lines[1].split(" ")

    hh, mm, ss = (int(x) for x in timestamp.split(":"))

    if year != int(_year):
        raise ValueError(
            f"Year mismatch in temperature log entry. Got {year} and {_year}"
        )
    if hour != hh:
        raise ValueError(f"Hour mismatch in temperature log entry. Got {hour} and {hh}")

    # don't care about the third line
    date_object = datetime(year, 1, 1, hh, mm, ss, tzinfo=timezone.utc) + timedelta(
        days=day - 1
    )
    try:
        return {
            "time": Time(date_object),
            "front_end_temperature": float(lines[3].split(" ")[1]) * un.deg_C,  # 100
            "amb_load_temperature": float(lines[4].split(" ")[1]) * un.deg_C,  # 101
            "hot_load_temperature": float(lines[5].split(" ")[1]) * un.deg_C,  # 102
            "inner_box_temperature": float(lines[6].split(" ")[1]) * un.deg_C,  # 103
            "thermal_control": float(lines[7].split(" ")[1]),  # 106
            "battery_voltage": float(lines[8].split(" ")[1]),  # 150
            "pr59_current": float(lines[9].split(" ")[1]),  # 152
        }
    except (ValueError, IndexError) as e:
        linestr = "".join(lines)
        raise ValueError(f"Error parsing temperature log entry:\n{linestr}") from e


def read_temperature_log(logfile: Path | str) -> QTable:
    """Read a full temperature log file."""
    pattern = re.compile(r"(\d{4})_(\d{3})_(\d{2})")
    with Path(logfile).open("r") as fl:
        lines = fl.readlines()

        # First, let's locate all the individual entries
        chunk_indices = []
        for i, line in enumerate(lines):
            if re.match(pattern, line) is not None:
                chunk_indices.append(i)
        chunk_indices.append(len(lines))

        # Only take chunks that are length 10
        chunk_indices = list(zip(chunk_indices[:-1], chunk_indices[1:]))
        chunk_indices = [c for c in chunk_indices if c[1] - c[0] == 10]

        all_data = []
        for c in chunk_indices:
            try:
                all_data.append(read_temperature_log_entry(lines[c[0] : c[1]]))
            except ValueError as e:
                warnings.warn(str(e), stacklevel=2)

    out = QTable(all_data)

    # Convert temperatures to Kelvin
    for col in out.columns:
        if col.endswith("temperature"):
            out[col] = out[col].to(un.K, equivalencies=un.temperature())

    return out


def get_mean_temperature(
    temperature_table: QTable,
    start_time: Time | None = None,
    end_time: Time | None = None,
    load: Literal["box", "amb", "hot", "open", "short"] = "box",
):
    """Get the mean temperature for a particular load from the temperature table."""
    if start_time is not None or end_time is not None:
        if start_time is None:
            start_time = temperature_table["time"].min()
        if end_time is None:
            end_time = temperature_table["time"].max()

        mask = (temperature_table["time"] >= start_time) & (
            temperature_table["time"] <= end_time
        )

        if not np.any(mask):
            raise ValueError(
                f"No data found between {start_time} and {end_time} in temperature table"
            )

        temperature_table = temperature_table[mask]

    if load in ("hot", "hot_load"):
        return temperature_table["hot_load_temperature"].mean()
    elif load in ("amb", "ambient", "open", "short", "box"):
        return temperature_table["amb_load_temperature"].mean()
    else:
        raise ValueError(f"Unknown load {load}")


@hickleable()
@attrs.define(frozen=True, hash=True)
class CalibrationObservation:
    """An EDGES-3 calibration observation."""

    s11_files: frozendict[str : frozendict[str, Path]] = attrs.field()
    acq_files: frozendict[str:Path] = attrs.field()
    temperature_file: Path = attrs.field()

    @classmethod
    def from_date(
        cls,
        year: int,
        day: int,
        s11_year: int | None = None,
        s11_day: int | None = None,
        s11_hour: int | str = "first",
        root_dir: Path | str = "/data5/edges/data/EDGES3_data/MRO/",
        allow_closest_s11_within: int = 5,
    ):
        """Create a CalibrationObservation from a date.

        Parameters
        ----------
        year : int
            The year of the observed spectra.
        day : int
            The day of the year of the observed spectra.
        s11_year : int, optional
            The year of the S11 files. If None, defaults to `year`.
        s11_day : int, optional
            The day of the S11 files. If None, defaults to `day`.
        s11_hour : int | str, optional
            The hour of the S11 files (if multiple files exist at different hours on
            the same day).
            If "first" or "last", will take the first or last file.
        root_dir : Path | str, optional
            The root directory of the EDGES-3 data.
        allow_closest_s11_within : int, optional
            The maximum number of days to search for the closest S11 file if none is
            found on the given day.
        """
        load_map = {
            "amb": "ambient",
            "hot": "hot_load",
            "open": "open",
            "short": "short",
            "lna": "lna",
        }
        if s11_year is None:
            s11_year = year
        if s11_day is None:
            s11_day = day

        s11_files: frozendict[str : frozendict[str, Path]] = frozendict(
            {
                load_map[load]: frozendict(
                    get_s1p_files(
                        load,
                        s11_year,
                        s11_day,
                        root_dir,
                        s11_hour,
                        allow_closest_s11_within,
                    )
                )
                for load in ["amb", "hot", "open", "short", "lna"]
            }
        )

        acq_files: dict[str:Path] = frozendict(
            {
                load_map[load]: get_acq_file(load, root_dir, year, day)
                for load in ["amb", "hot", "open", "short"]
            }
        )

        temperature_file = root_dir / "temperature_logger/temperature.log"

        return cls(s11_files, acq_files, temperature_file)

    def get_spectra(self, load: str) -> io.FieldSpectrum:
        return io.FieldSpectrum(self.acq_files[load])

    @cache  # noqa: B019
    def get_temperature_table(self):
        return read_temperature_log(self.temperature_file)

    def __gethstate__(self):
        return attrs.asdict(self)
