from __future__ import annotations

import datetime
import numpy as np
import os
import shutil
from itertools import islice
from pathlib import Path
from rich.console import Console
from typing import List, Optional, Union

IGNORABLE = (".old", ".ignore", ".invalid", ".output")
console = Console()


def make_symlink_tree(files: dict[str, Path], symdir: Path, obs_name):

    for fl, fl_abs in files.items():
        sym_path = symdir / obs_name / fl

        # Here we make all the directories, that all have to be symlinks, otherwise
        # objects whose path is to a directory don't equate.
        if not sym_path.parent.exists():
            sym_path.parent.mkdir(parents=True)

        sym_path.symlink_to(fl_abs)


def get_active_files(path: str | Path) -> list[Path]:
    path = Path(path)
    if not path.is_dir():
        raise ValueError(f"{path} is not a directory!")
    fls = path.glob("*")
    ok_extra_files = [
        "Notes.txt",
        "calibration_analysis.ipynb",
        "derived",
        "definition.yaml",
    ]

    return [
        fl for fl in fls if fl.suffix not in IGNORABLE and fl.name not in ok_extra_files
    ]


def get_parent_dir(path, n=1):
    for _ in range(n):
        path = os.path.dirname(os.path.normpath(path))
    return path


def ymd_to_jd(y, m, d):
    return (
        datetime.date(int(y), int(m), int(d)) - datetime.date(int(y), 1, 1)
    ).days + 1


def _ask_to_rm(fl: Path) -> Path | None:
    reply = "z"
    while reply not in "yniopm":
        reply = (
            str(
                console.input(
                    f"""[bold] Action for '{fl}':
>>> [red]y[/]: remove the file
>>> [dark_orange]n[/]: ignore/skip the file for now
>>> [cyan]i[/]: make it 'invalid'
>>> [cyan]o[/]: [bold]make it 'old'
>>> [cyan]p[/]: [bold]make it 'output'
>>> [green]m[/]: [bold]interactively move/rename
"""
                )
            )
            .lower()
            .strip()
        )

    if reply == "y":
        if fl.is_dir():
            shutil.rmtree(fl)
        else:
            os.remove(str(fl))
        return None
    elif reply == "i":
        shutil.move(fl, str(fl) + ".invalid")
        return Path(str(fl) + ".invalid")
    elif reply == "o":
        shutil.move(fl, str(fl) + ".old")
        return Path(str(fl) + ".old")
    elif reply == "p":
        shutil.move(fl, str(fl) + ".output")
        return Path(str(fl) + ".output")
    elif reply == "m":
        reply = str(console.input(f"[bold]Change '{fl.name}' to:[/] "))
        newfile = fl.parent / reply
        try:
            shutil.move(fl, newfile)
        except Exception:
            console.print(
                f"[bold red]Couldn't rename the file '{newfile}' as you asked."
            )
            raise
        return newfile
    elif reply == "n":
        return fl
    else:
        raise OSError("Something went very wrong.")


class FileStructureError(Exception):
    pass


class LoadExistError(Exception):
    pass


class IncompleteObservation(FileStructureError):
    pass


class InconsistentObservation(FileStructureError):
    pass


def get_file_list(top_level: Path, filter=None, ignore=None):
    ignore = ignore or []

    out = []
    for pth in top_level.iterdir():
        if str(pth) not in ignore and (filter(pth) if filter is not None else True):
            if pth.is_file():
                out.append(pth.absolute())
            elif pth.is_dir():
                out.extend(get_file_list(pth, filter=filter, ignore=ignore))
    return out


def snake_to_camel(word: str):
    return "".join(w[0].upper() + w[1:] for w in word.split("_"))


def optional(fnc):
    def outfnc(x):
        return x is None or fnc(x)

    outfnc.optional = True
    return outfnc


def isintish(x):
    return isinstance(x, (int, np.int64))


def isfloatish(x):
    return isinstance(x, (float, np.float32))


def isstringish(x):
    return isinstance(x, (str, bytes))


def isnumeric(x):
    return isinstance(x, (int, float, np.number))
