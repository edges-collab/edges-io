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

As of v0.4.0, ``CalibrationObservation`` objects no longer need to be defined fully by
one directory containing all measurements. While that is still an option (and the easiest
way to define a calibration observation), they can also be defined in a more sophisticated
way internally or externally.

Internally, a ``definition.yaml`` file is allowed which defines properties of the
observation, and also has ``include`` and ``prefer`` keywords which are used to supplement
or override any particular parts of the observation. For example ``include`` could point
to the top-level of any other observation, which could then be used whenever the
main observation lacks data. If this file exists, by default it is used to construct
the full observation virtually. See an incomplete example of such a definition file can
be found `here <example-obs-definition.yaml>`_.


Externally, a different file format is used to explicitly define every single measurement
file in an observation. This is supposed to be exhaustive and complete to make it
unambiguous. An example can be found in the `test-suite <tests/test_data/observation.yaml>`_.
One can use such a file to create a ``CalibrationObservation`` by using the
``CalibrationObservation.from_observation_yaml()`` function.

Using the ``HDF5Object``
------------------------
``edges-io`` contains a convenient ``HDF5Object`` class whose purpose is to make working
with HDF5 files a bit more formal (and arguably more simple). By subclassing it, you
can specify an exact file layout that can be verified on read, to ensure a file is
in the correct format (not just HDF5, but that it has the correct data structures and
groups and metadata).

Using such a class is meant to provide a very thin wrapper over the file. So, for instance
if you have a file ``my_hdf5_format_file.h5``, whose structure is defined by the class
``CustomH5Format``, you can create an object like this::

    >>> fl = CustomH5Format("my_hdf5_format_file.h5")

Directly on creation, the file will be checked for compatibility and return an error
if it contains extraneous keys, or lacks keys that it requires.

Once created, the ``fl`` variable now has operations which can "look into" the file
and load its data. It supports lazy-loading, so doing::

    >>> print(fl['dataset'].max())

will load the 'dataset' data, and get the maximum, but it will not keep the data in
memory, and will not load any other datasets. If you have data in groups, you can
easily do::

    >>> print(fl['group']['dataset'].min())

To load the data into the object permanently use the ``.load`` method::

    >>> fl.load('group')

In fact, doing this will load all data under 'group'. If you just wanted to load
"dataset" out of "group"::

    >>> fl['group'].load('dataset')

An example of how to define a subclass of ``HDF5Object`` can be seen in the
``HDF5RawSpectrum`` class, which is used to define fastspec output files.


Note
====

This project has been set up using PyScaffold 3.2.3. For details and usage
information on PyScaffold see https://pyscaffold.org/.
