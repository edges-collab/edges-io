import re
import shutil
import subprocess
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Iterable, Tuple, Union

from . import utils
from .logging import logger


class _ObsNode(ABC):
    """Abstract base class representing a node in a calibration observation.

    A node could be a file or a directory.
    """

    known_substitutions = ()
    known_patterns = ()
    pattern = ""
    write_pattern = ""

    def __init__(
        self,
        path: [str, Path],
        *,
        check: bool = True,
        fix: bool = False,
        log_level: int = 40,
    ):

        if check:
            pre_level = logger.getEffectiveLevel()
            logger.setLevel(log_level)

            self.path, self._match_dict = self._run_checks(path, fix)

            logger.setLevel(pre_level)
        else:
            self.path = Path(path)
            self._match_dict = None

    @classmethod
    def _run_checks(cls, path, fix):
        return cls.check_self(path, fix)

    @staticmethod
    @abstractmethod
    def check_self(
        path, fix=False
    ) -> Tuple[Path, Union[dict[str, str], Iterable[dict[str, str]]]]:
        pass

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
    def _fix(cls, root: Path, basename: str):
        """Auto-fix a basename."""

        # First try simple substitutions
        new_name = basename
        for sub, correct in cls.known_substitutions:
            if sub in new_name:
                new_name = new_name.replace(sub, correct)

        match = re.search(cls.pattern, new_name)

        # If a simple substitution did the trick, return.
        if match is not None:
            shutil.move(root / basename, root / new_name)
            logger.success(f"Successfully converted to {new_name}")
            return root / new_name, match

        # Otherwise, try various patterns.
        for pattern in cls.known_patterns:
            match = re.search(pattern, new_name)
            if match:
                break

        if match is None:
            logger.warning("\tCould not auto-fix it.")
            fixed = utils._ask_to_rm(root / basename)
            if fixed:
                logger.success("Successfully removed.")
            return None, None
        else:
            dct = match.groupdict()
            default = {
                **cls._get_filename_parameters(dct),
                **cls._get_filename_params_from_contents(root / new_name),
            }
            dct = {**default, **dct}

            new_name = cls.write_pattern.format(**dct)
            new_path = root / new_name

            match = re.search(cls.pattern, new_name)

            if match is not None:
                logger.success(f"Successfully converted to {new_name}")
                shutil.move(root / basename, new_path)
                return new_path, match
            else:
                return None, None


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


class _DataContainer(_ObsNode):
    _content_type = None

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
    def check_contents(cls, path: [str, Path], fix=False) -> bool:
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
    def _check_contents_selves(cls, path: Path, fix=False) -> bool:
        fls = utils.get_active_files(path)

        # Start off with a clean slate for this function.
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
                                    [f"pandoc -o {fl.with_suffix('.txt')} {fl}"],
                                    check=True,
                                )
                                fl = fl.with_suffix(".txt")
                                logger.success(f"Successfully renamed to {fl}")
                            except subprocess.CalledProcessError as e:
                                logger.warning(
                                    f"Could not convert to .txt -- error: {e.message}"
                                )

                        else:
                            fixed = utils._ask_to_rm(fl)

                            if fixed:
                                logger.success("Successfully removed.")

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
        logger.errored = 0
        return ok

    @classmethod
    def typestr(cls, name: str) -> str:
        """Generate a string uniquely defining the 'kind of thing' path is.

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
    def _fix(cls, root, basename):
        return
