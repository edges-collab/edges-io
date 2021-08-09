=========
Changelog
=========
v2.5.0
======

Added
-----

* Ability to register types as "YAML", so that they can be read from the HDF5 ``attrs``
  (or "meta") as native types via YAML. This works both ways (i.e. read and write) as
  long as a constructor/representer is defined for the object in YAML.

v2.4.0
======

Added
-----
* Ability to open ``HDF5Object`` instances' files as writeable, so as to append extra data
  to them.

v2.3.9
======

Try again to get this on PyPI... grrr :-(

v2.3.7
======

A whole bunch of minor updates just to get the auto-deployment working...

v2.3.1
======

First version on PyPI!

v2.0.0
======
Changed
-------
* Removed ``.load`` and ``.load_all`` methods from ``HDF5Object`` so that the interface
  between the object and its sub-groups are more identical. Note that *everything* now
  caches when accessed, but can be explicitly removed from the cache at will.

v0.5.0
======
Added
-----
* ``list_of_files`` method that contains all files used in a calibration.

v0.4.0
======
Added
-----
* ``observation.yaml`` file for explicitly defining a full calibration observation.
* ``definition.yaml`` file for internally defining metadata of an observation, with the
  ability to include other observations to supplement/override it.

Changed
-------
* Default behaviour when instantiating ``CalibrationObservation`` is now to not print
  out structure and info logs, only errors.

v0.3.0
======
Added
-----
* TOML file for antenna simulators
* Access to all Resistance metadata via a structured numpy array ``raw_data``.
* Antenna simulators are now able to be identified and read in more easily from each component.

Changed
-------
* Better way of getting the version
* Ignore some more different kinds of files.

Fixed
-----
* New checks on whether antenna simulators are the same for S11, Spectra and Resistance.

v0.2.0
======

Added
-----
* New auto-fixes for the root directory and spectra.

Changed
-------
* load_name is now a simpler alias (ambient, hot_load, open, short)

Fixed
-----
* Several bug-fixes found when fixing actual data on disk.


v0.1.0
======

- First public version
