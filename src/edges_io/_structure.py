from __future__ import annotations

import attr
import os
import re
import shutil
import subprocess
from abc import ABC, abstractmethod
from cached_property import cached_property
from copy import copy
from pathlib import Path
from typing import Dict, Iterable, Optional, Tuple, Union

from . import utils
from .logging import logger


@attr.s(frozen=True)
class _ObsNode(ABC):
    """Abstract base class representing a node in a calibration observation.

    A node could be a file or a directory.
    """

    known_substitutions: tuple[str] = ()
    known_patterns: tuple[str] = ()
    pattern: str = ""
    write_pattern: str = ""

    _path: str | Path = attr.ib(converter=Path)
    check: bool = attr.ib(default=True, kw_only=True, converter=bool)
    fix: bool = attr.ib(default=False, kw_only=True, converter=bool)
    _log_level: int = attr.ib(default=40, converter=int)

    def __attrs_post_init__(self):
        # Here we just access this to ensure the checks run on instantiation.
        self._path_and_match_dict

    @cached_property
    def _path_and_match_dict(self) -> tuple[Path, dict | None]:
        if self.check:
            pre_level = logger.getEffectiveLevel()
            logger.setLevel(self._log_level)
            path, match_dict = self._run_checks(self._path, self.fix)
            logger.setLevel(pre_level)
            return path, match_dict
        else:
            return self._path, None

    @property
    def path(self) -> Path:
        return self._path_and_match_dict[0]

    @property
    def _match_dict(self) -> dict | None:
        return self._path_and_match_dict[1]

    @classmethod
    def _run_checks(cls, path, fix):
        return cls.check_self(path, fix)

    @classmethod
    def check_self(cls, path: str | Path, fix: bool) -> tuple[Path, dict | None]:
        path = Path(path)
        path, match = cls._check_self(path, fix=fix)
        if match is not None:
            cls._validate_match(match, path.name)
        return path, match

    @classmethod
    def _check_self(cls, path: Path, *, fix: bool = False) -> tuple[Path, dict | None]:

        if not path.exists():
            logger.error(f"The path {path} does not exist!")

        match = re.search(cls.pattern, path.name)
        if match is None:
            logger.error(
                f"The filename {path.name} does not have the correct format for a {cls.__name__}. "
                f"Correct format: {cls.write_pattern}."
            )

            if fix:
                path, match = cls._fix(path.parent, path.name)

        return path, match.groupdict() if match is not None else None

    @classmethod
    def typestr(cls, name: str):
        """Generate a string uniquely defining the 'kind of thing' the object is.

        The point of this method is to be able to compare two different file/folder
        names to check whether they describe the same kind of thing. For example,
        two Spectrum files from different observations which are both "Ambient" should
        return the same string, even though their dates etc. might be different. However,
        two Spectrum files of different Loads will be different.

        The reason this has to exist (and as a classmethod) is because merely comparing
        the type of the thing is not enough, since the class itself has no knowledge of
        for example the kind of load. But comparing instances is not great either, since
        instances are expected to have a full complement of files to be "valid", but
        one of the main purposes of comparing files is to construct such a full observation.
        """
        return cls.__name__

    @classmethod
    def _get_filename_parameters(cls, dct: dict):
        """Return a dictionary of filename parameters to be inserted into `write_pattern`.

        These should be defaults if appropriate.
        If a particular parameter has no default and cannot be obtained, omit it."""
        return {}

    @classmethod
    def _get_filename_params_from_contents(cls, path: Path):
        return {}

    @classmethod
    def _fix(cls, root: Path, basename: str) -> tuple[Path | None, re.Match | None]:
        """Auto-fix a basename."""

        # First try simple substitutions
        new_name = copy(basename)
        for sub, correct in cls.known_substitutions:
            if sub in new_name:
                new_name = new_name.replace(sub, correct)

        match = re.search(cls.pattern, new_name)

        # If a simple substitution did the trick, return.
        if new_name != basename:
            shutil.move(root / basename, root / new_name)

        if match is not None:
            logger.success(f"Successfully converted to {new_name}")
            return root / new_name, match

        basename = copy(new_name)

        # Otherwise, try various patterns.
        for pattern in cls.known_patterns:
            match = re.search(pattern, new_name)
            if match:
                break

        if match is None:
            logger.warning(f"Could not auto-fix {basename}.")
            new_path = utils._ask_to_rm(root / basename)

            if new_path is None:
                logger.success("Successfully removed.")
            elif new_path != (root / basename):
                if new_path.suffix not in utils.IGNORABLE:
                    match = re.search(cls.pattern, str(new_path.relative_to(root)))

                    if match is None:
                        logger.warning(
                            f"New name '{new_path.relative_to(root)}' is not correctly "
                            "formatted either!"
                        )

                    else:
                        logger.success(f"Successfully moved to '{new_path.name}'")
                else:
                    logger.success(f"Successfully moved to '{new_path.name}'")

                basename = new_path.name

            return root / basename, match
        else:
            dct = match.groupdict()
            dct = {
                **cls._get_filename_parameters(dct),
                **dct,
                **cls._get_filename_params_from_contents(root / basename),
            }

            new_name = cls.write_pattern.format(**dct)
            new_path = root / new_name

            match = re.search(cls.pattern, new_name)

            if match is None:
                return root / basename, None
            logger.success(f"\tSuccessfully converted to {new_name}")
            shutil.move(root / basename, new_path)
            return new_path, match

    @classmethod
    def _validate_match(cls, match: dict[str, str], filename: str):
        pass


@attr.s
class _DataFile(_ObsNode):
    """
    Abstract Object representing a file in a calibration observation.

    Parameters
    ----------
    path
        The path to the file.
    fix
        Whether to attempt to fix the file in place if its filename is in the wrong
        format.
    """

    @classmethod
    def _check_self(cls, path: Path, **kwargs) -> tuple[Path, dict | None]:
        return super()._check_self(path, **kwargs)


@attr.s
class _DataContainer(_ObsNode):
    _content_type = None

    @classmethod
    def _check_self(cls, path: Path, **kwargs) -> tuple[Path, dict | None]:
        logger._structure(f"Checking {cls.__name__} folder contents format at {path}.")
        return super()._check_self(path, **kwargs)

    @classmethod
    def _run_checks(cls, path, fix):

        path, match = cls.check_self(path, fix)

        if match is None:
            raise utils.FileStructureError(
                f"Directory {path.name} is in the wrong format."
            )

        if not cls._check_contents_selves(path, fix):
            raise utils.FileStructureError()
        if not cls._check_all_files_there(path):
            raise utils.IncompleteObservation()
        if not cls._check_file_consistency(path):
            raise utils.InconsistentObservation()

        return path, match

    @classmethod
    def check_contents(cls, path: str | Path, fix=False) -> bool:
        """Abstract method for checking whether the contents of this container are in
        the correct format for the DB"""
        # Check that everything that *is* there has correct format.
        path = Path(path)
        ok_selves = cls._check_contents_selves(path, fix=fix)
        ok_complete = cls._check_all_files_there(path)
        # Check that the files that are there have consistent properties, and are also
        # consistent with outside parameters (eg. if year appears on them, they should
        # be consistent with outer years).
        ok_consistent = cls._check_file_consistency(path)

        return ok_selves and ok_complete and ok_consistent

    @classmethod
    @abstractmethod
    def _check_all_files_there(cls, path: Path) -> bool:
        return True

    @classmethod
    @abstractmethod
    def _check_file_consistency(cls, path: Path) -> bool:
        return True

    @classmethod
    def _check_contents_selves(cls, path: Path, fix=False) -> tuple[bool, Path]:
        fls = utils.get_active_files(path)

        # Start off with a clean slate for this function.
        n_errors = logger.errored
        logger.errored = 0
        for fl in fls:
            if isinstance(cls._content_type, dict):
                for key, ct in cls._content_type.items():
                    if fl.name.startswith(key):
                        content_type = ct
                        break
                else:
                    logger.error(f"{fl.name} is an extraneous file/folder")

                    if fix:
                        if fl.name == "Notes.odt":
                            try:
                                subprocess.run(
                                    [
                                        "pandoc",
                                        "-o",
                                        str(fl.with_suffix(".txt")),
                                        str(fl),
                                    ],
                                    check=True,
                                )
                                fl = fl.with_suffix(".txt")
                                if fl.exists():
                                    os.remove(str(fl.with_suffix(".odt")))

                                logger.success(f"Successfully converted to {fl}")
                            except subprocess.CalledProcessError as e:
                                logger.warning(
                                    f"Could not convert to .txt -- error: {e.message}"
                                )

                        else:
                            new_path = utils._ask_to_rm(fl)

                            if new_path is None:
                                logger.success("Successfully removed.")
                            elif new_path != fl:
                                if new_path.suffix not in utils.IGNORABLE:
                                    match = re.search(
                                        cls.pattern,
                                        str(new_path.relative_to(fl.parent)),
                                    )

                                    if match is None:
                                        logger.warning(
                                            f"New name '{str(new_path.relative_to(fl.parent))}' is not "
                                            f"correctly formatted either!"
                                        )
                                    else:
                                        logger.success(
                                            f"Successfully moved to '{new_path.name}'"
                                        )
                                else:
                                    logger.success(
                                        f"Successfully moved to '{new_path.name}'"
                                    )

                    continue
            else:
                content_type = cls._content_type

            fl, _ = content_type.check_self(fl, fix=fix)

            # Recursively check the contents of the contents.
            try:
                content_type.check_contents(fl, fix=fix)
            except AttributeError:
                # It's a DataFile, not a DataContainer
                pass

        ok = not bool(logger.errored)
        logger.errored = n_errors + logger.errored
        return ok
