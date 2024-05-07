"""Simple type definitions for use internally."""

from pathlib import Path
from typing import Union

from astropy import units

PathLike = Union[str, Path]
FreqType = units.Quantity["frequency"]
ImpedanceType = units.Quantity["electrical impedance"]
OhmType = units.Quantity[units.ohm]
TimeType = units.Quantity["time"]
DimlessType = units.Quantity["dimensionless"]
