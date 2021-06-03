from pkg_resources import DistributionNotFound, get_distribution

try:
    __version__ = get_distribution(__name__).version
except DistributionNotFound:
    __version__ = "unknown"
finally:
    del get_distribution, DistributionNotFound

from . import h5, io
from .h5 import HDF5Object, HDF5RawSpectrum, HDF5StructureError
from .io import S1P, FieldSpectrum, Spectrum
