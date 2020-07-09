=========
Changelog
=========

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
