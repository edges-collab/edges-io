import datetime
import os
import shutil
from pathlib import Path
from typing import List


def get_active_files(path: [str, Path]) -> List[Path]:
    path = Path(path)
    if not path.is_dir():
        raise ValueError(f"{path} is not a directory!")
    fls = path.glob("*")
    ok_extra_files = ["Notes.txt", "calibration_analysis.ipynb", "derived"]

    return [
        fl
        for fl in fls
        if fl.suffix not in (".old", ".ignore", ".invalid", ".output")
        and fl.name not in ok_extra_files
    ]


def get_parent_dir(path, n=1):
    for _ in range(n):
        path = os.path.dirname(os.path.normpath(path))
    return path


def ymd_to_jd(y, m, d):
    return (
        datetime.date(int(y), int(m), int(d)) - datetime.date(int(y), 1, 1)
    ).days + 1


def _ask_to_rm(fl: Path):
    while True:
        reply = (
            str(
                input(
                    f"Would you like to (recursively) remove {fl} ([y]es/make [i]nvalid/make [o]ld/make out[p]ut/[m]ove/[n]o)?: "
                )
            )
            .lower()
            .strip()
        )
        if reply.startswith("y"):
            rm = True
            break
        elif reply.startswith("n") or not reply:
            rm = False
            break
        elif reply.startswith("i"):
            rm = None
            kind = "i"
            break
        elif reply.startswith("o"):
            rm = None
            kind = "o"
            break
        elif reply.startswith("p"):
            rm = None
            kind = "p"
            break
        elif reply.startswith("m"):
            rm = None
            kind = "m"
            break
        else:
            print("please select (y/n) only")

    if rm:
        if fl.is_dir():
            shutil.rmtree(fl)
        else:
            os.remove(str(fl))
        return True
    elif rm is None:
        if kind == "i":
            shutil.move(fl, str(fl) + ".invalid")
            return True
        if kind == "o":
            shutil.move(fl, str(fl) + ".old")
            return True
        if kind == "p":
            shutil.move(fl, str(fl) + ".output")
            return True
        elif kind == "m":
            reply = str(input(f"Change {fl.name} to: "))
            newfile = fl.parent / reply
            try:
                shutil.move(fl, newfile)
            except Exception:
                print(f"Couldn't rename the file {newfile} as you asked.")
                raise
    else:
        return False


class FileStructureError(Exception):
    pass


class IncompleteObservation(FileStructureError):
    pass


class InconsistentObservation(FileStructureError):
    pass


def get_file_list(top_level: Path, filter=None, ignore=None):
    ignore = ignore or []

    out = []
    for pth in top_level.iterdir():
        if (
            pth.is_file()
            and str(pth) not in ignore
            and (filter(pth) if filter is not None else True)
        ):
            out.append(pth.absolute())
        elif pth.is_dir():
            out.extend(get_file_list(pth, filter=filter, ignore=ignore))
    return out


def snake_to_camel(word: str):
    return "".join(w.capitalize() for w in word.split("_"))
