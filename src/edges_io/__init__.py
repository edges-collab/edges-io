from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version(__name__)
except PackageNotFoundError:
    # package is not installed
    __version__ = "unknown"
finally:
    del version, PackageNotFoundError

from pathlib import Path

from . import io, io3
from .io import S1P, FieldSpectrum, Spectrum

TEST_DATA_PATH = Path(__file__).parent / "test_data"

del Path
