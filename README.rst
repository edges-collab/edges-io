========
edges-io
========

.. image:: https://travis-ci.org/edges-collab/edges-io.svg?branch=master
    :target: https://travis-ci.org/edges-collab/edges-io
.. image:: https://codecov.io/gh/edges-collab/edges-io/branch/master/graph/badge.svg
  :target: https://codecov.io/gh/edges-collab/edges-io

**Module for reading EDGES data and working with EDGES databases.**

This package implements all necessary functionality for reading EDGES data.
It's two main concerns are:

1. Reading the various file formats required for VNA readings, fastspec, etc.
2. Verifying and exploring databases of measurements in a robust and reliable way.

Features
========
Some features currently implemented:

* Verify a "calibration observation" quickly without reading any actual data, with
  a nice command-line tool: ``edges-io check``.
* Optionally apply various automatic _fixes_ to a calibration observation to bring
  it into line with standard database layout.
* Read ``acq``, ``h5``, ``mat`` and ``npz`` spectrum files seamlessly.
* Read S1P files.
* Verification of read data.
* Intuitive class hierarchy so that any subset of an observation can be handled.

Installation
============
Installation should be as simple as either one of the following::

    $ pip install git+git://github.com/edges-collab/edges-io

or, if you would like to develop ``edges-io`` and use it too::

    $ git clone https://github.com/edges-collab/edges-io
    $ cd edges-io
    $ pip install -e .[dev]

There are a few dependencies, which should be installed automatically when following the
above command. If you are using ``conda`` (which is recommended) then you can obtain
a cleaner/faster install by doing the following::

    $ conda create -n edges python=3
    $ conda activate edges
    $ conda install numpy scipy h5py

And then following either of the above instructions.

Usage
=====
You can use ``edges-io`` either as a library or a command-line tool. The library is
self-documented, so you can look at the docstring of any of the available functions.
We describe some basics of each approach here.

CLI
---
To run the checking tool, simply do::

    $ edges-io check PATH

``PATH`` should be the top-level directory of a calibration observation (i.e. a folder
that has a sub-folder ``25C/``, which has subfolders ``Spectra/``, ``Resistance/`` and
``S11`` etc.).
There are a few options you can use, example changing the temperature of the observation,
and enabling automatic fixes. The latter can be achieved simply with the ``--fix`` flag.
If you find that a particular kind of error happens regularly,
`make an issue <https://github.com/edges-collab/edges-io/issues/new>`_ so we can add the
fix.

Library
-------
The library is useful for gathering an entire observation and performing operations
on its data. The library exposes a hierarchy of calibration objects, including base
objects like a ``Spectrum``, ``Resistance`` or ``S1P`` file, and container objects
like ``Spectra`` or ``S11``. An entire observation can be loaded as a
``CalibrationObservation``, and it contains references to all children.

For example::

    >>> from edges_io import io
    >>> obs = io.CalibrationObservation("path_to_observation")
    >>> obs.s11.path
    "path_to_observation/25C/S11"
    >>> obs.spectra.ambient.path
    "path_to_observation/25C/Spectra/Ambient_XXX.acq"
    >>> ambient_spectrum = obs.spectra.ambient.read()

See how ``edges-io`` is used in
`edges-cal <https://github.com/edges-collab/cal_coefficients/tree/master/src/edges_cal/cal_coefficients.py>`_
for a more involved example.

Note
====

This project has been set up using PyScaffold 3.2.3. For details and usage
information on PyScaffold see https://pyscaffold.org/.
