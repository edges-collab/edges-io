import datetime
import glob
import os
import shutil
from pathlib import Path


def get_active_files(path):
    if not os.path.isdir(path):
        raise ValueError("{} is not a directory!".format(path))
    fls = glob.glob(os.path.join(path, "*"))
    ok_extra_files = ["Notes.txt", "calibration_analysis.ipynb", "derived"]

    return [
        fl
        for fl in fls
        if not fl.endswith(".old")
        and os.path.basename(fl) != "Notes.txt"
        and not fl.endswith(".ignore")
        and not fl.endswith(".invalid")
        and os.path.basename(fl) not in ok_extra_files
    ]


def get_parent_dir(path, n=1):
    for _ in range(n):
        path = os.path.dirname(os.path.normpath(path))
    return path


def ymd_to_jd(y, m, d):
    return (
        datetime.date(int(y), int(m), int(d)) - datetime.date(int(y), 1, 1)
    ).days + 1


def _ask_to_rm(fl):
    while True:
        reply = (
            str(
                input(
                    "Would you like to (recursively) remove {} ([y]es/[i]gnore/[m]ove/[n]o)?: ".format(
                        fl
                    )
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
        elif reply.startswith("m"):
            rm = None
            kind = "m"
            break
        else:
            print("please select (y/n) only")

    if rm:
        if os.path.isdir(fl):
            shutil.rmtree(fl)
        else:
            os.remove(fl)
        return True
    elif rm is None:
        if kind == "i":
            shutil.move(fl, fl + ".old")
            return True
        elif kind == "m":
            reply = str(input("Change {} to: ".format(os.path.basename(fl))))
            newfile = os.path.join(os.path.dirname(os.path.normpath(fl)), reply)
            try:
                shutil.move(fl, newfile)
            except Exception:
                print("Couldn't rename the file {} as you asked.".format(newfile))
                raise
    else:
        return False


class FileStructureError(Exception):
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
