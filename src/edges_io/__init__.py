# -*- coding: utf-8 -*-
from pkg_resources import get_distribution, DistributionNotFound

try:
    __version__ = get_distribution(__name__).version
except DistributionNotFound:
    __version__ = "unknown"
finally:
    del get_distribution, DistributionNotFound

from . import io
from .io import Spectrum, S1P
from . import h5
from .h5 import HDF5Object, HDF5RawSpectrum, HDF5StructureError
