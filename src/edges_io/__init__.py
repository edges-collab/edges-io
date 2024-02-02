from pkg_resources import DistributionNotFound, get_distribution

try:
    __version__ = get_distribution(__name__).version
except DistributionNotFound:
    __version__ = "unknown"
finally:
    del get_distribution, DistributionNotFound

from pathlib import Path

from . import h5, io, io3
from .h5 import HDF5Object, HDF5RawSpectrum, HDF5StructureError
from .io import S1P, FieldSpectrum, Spectrum

TEST_DATA_PATH = Path(__file__).parent / "test_data"

del Path
