========
edges-io
========

.. image:: https://travis-ci.org/edges-collab/edges-io.svg?branch=master
    :target: https://travis-ci.org/edges-collab/edges-io
.. image:: https://codecov.io/gh/edges-collab/edges-io/branch/master/graph/badge.svg
  :target: https://codecov.io/gh/edges-collab/edges-io
.. image:: https://readthedocs.org/projects/edges-io/badge/?version=latest
  :target: https://edges-io.readthedocs.io/en/latest/?badge=latest
  :alt: Documentation Status

**Module for reading EDGES data and working with EDGES databases.**

This package implements all necessary functionality for reading EDGES data.
It's two main concerns are:

1. Reading the various file formats required for EDGES data:
   - VNA readings
   - fastspec output
   - thermistor readings
   - field weather recordings
   - field thermlog recordings
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
* Read field-based weather and thermlog information

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
There are a few options you can use, for example changing the temperature of the observation,
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
    >>> print(obs.s11.path)
    "path_to_observation/25C/S11"
    >>> print(obs.spectra.ambient.path)
    "path_to_observation/25C/Spectra/Ambient_XXX.acq"
    >>> ambient_spectrum = obs.spectra.ambient.read()

See how ``edges-io`` is used in
`edges-cal <https://github.com/edges-collab/cal_coefficients/tree/master/src/edges_cal/cal_coefficients.py>`_
for a more involved example.

Defining Observations
---------------------
One of the main goals of ``edges-io`` is to make the definition of a "Calibration
Observation" as clear, robust and error-free as possible. Many files go into any
particular observation -- spectra, resistance measurements and S11 measurements -- which
are all required to form together a calibration solution (which can then be applied
independently to field data). This code provides a clear structure for how these files
*must* be laid out in order for them to be read and used automatically. This is done in
a formal sense in `this document <docs/structure.rst>`_, but is also implemented within
the code itself.

In the above document, the specification is laid out as formally as possible, and that
document has the final word on what is allowable. However, this can mean it's a bit hard
to interpret, and so we here present a "simpler" guide to what constitutes a "Calibration
Observation".

A single natural "Observation" (see below for how to combine multiple observations into
a single "virtual" observation) is a single directory with multiple files/subdirectories
in it. That directory must be named under a certain convention that time-stamps it and
gives some useful metadata of the observation (like which receiver number was measured).
It is possible that in the future, a metadata file within the observation will specify
most of this information, but it is also useful to have a unique label for the observation.

One question that is important in all of this definition is what to do when either 1)
a file exists that shouldn't be there, or 2) a file doesn't exist that *should* be
there. It may be tempting to overlook extra files that shouldn't be there. However,
they can be a source of error. For example, spectra can be split across multiple files,
and we use a file pattern to find the files that should be read in. If an extra file
that "shouldn't be there" exists and the file pattern matches it, then errors can occur
(even worse if the contents of that file are able to be read by the spectrum reader, but
correspond to a different load or something of that sort, where the results will be wrong,
but no error raised). Thus, when *checking* the integrity of an observation we flag
extra files as *errors*, and require the user to fix them up. To make this a bit easier,
and let those files stay in the directory (so we don't lose potentially valuable information),
one of a few extensions can be added to the extra file:

  * ``.old``: for files that contain valid data but that is superseded by newer measurements
    and should be ignored,
  * ``.invalid``: for data that has something wrong with it (equipment broken, wrong input
    parameters, etc.),
  * ``.ignore``: files to ignore for any other reason.

If the file does not have one of these extensions, and is not in the list of accepted
files for the Observation, an error will be raised by the checker.

On the other hand, if a file is *missing* that must be there, different things can happen
in different situations. The default case is to treat this as a warning, which may be
counter-intuitive (surely missing a required file should be an error?!). The reason for
this is that that file may be supplemented by a different Observation. Perhaps this
Observation is incomplete  -- maybe all the data that was taken was a single set of
Spectra, which is supposed to complement a previous observation which had a full set
of measurements. In this case, while the "natural" Observation is incomplete, it is not
necessary to give an error, as long as a warning is given such that it must be combined
with another observation. Nevertheless, some combinations of files are required to have
been taken in the same physical observation to ensure consistency (namely, S11 measurements
for each standard in a given load). If particular standards are missing, an error will
be raised.

These caveats should be kept in mind as we talk about "required" directories/files below.
"Required" will mean that after combining all the observations that we want/need (see
next section), we require this particular file.

Within the top-level observation directory are a number of directories denoting the ambient temperature at
which the observation was taken. These will usually be 15C, 25C or 35C. Most newer
observations are at 25C. One should never mix files between different ambient temperatures.
Thus, in reality, an observation is contained within one of *these* folders, and in practice,
the ``CalibrationObservation`` has its ``path`` attribute set to the temperature directory.

Inside this directory can be up to two files, and exactly three folders. One of the files
is a ``Notes.txt`` file which summarises human-readable notes about the observation ("we
ran the ambient spectra first, but had a delay because of xxx..."). The other file is named
``definition.yaml`` and includes metadata about the observation in a specific format
(this file also allows you to supplement the observation with other observations, but
we'll get to that later). Measurements/data like the male/female resistance should be
put in here (til now they have been found somewhere an input manually by the analyst
when doing the calibration, which is very risky and prone to error -- they are properly
part of the measurement itself, not a choice of the analyst).

The three folders are ``Spectra``, ``Resistance`` and ``S11``. Note that an observation
must have all three of these (and nothing else, after combining observations).

Within ``Spectra`` exist a bunch of spectra taken over about 12-24 hours for each of
four "calibration sources" in the lab: they are "Ambient", "HotLoad", "LongCableOpen"
and "LongCableShorted" (often referred to as their simple aliases "ambient", "hot_load",
"open" and "short" in the code). These spectra will be in either ``.acq`` or ``.h5``
format, depending on the version of ``fastspec`` that took the measurements. Due to the
way ``fastspec`` takes its data, each source may have multiple files for a single
measurement (each integration is saved to a new line in the file, but a new file is
created at particular local times each day). Thus, typically one would like to read in
and concatenate _all_ the files for that load, to use all the data.
Beyond this, it is _possible_ that two fully separate "runs" for a given source/load
will be made. In this case, an identifier for the "run number" is put into the filename.
Only one run number is actually used to do any particular calibration. In practice, it
is *very rare* to have more than one "useable" run number for any particular load.
Typically, a second run is only taken if it is deemed necessary due to the first being
invalid in some way. If this is the case, this should be noted in the ``Notes.txt``
and/or the ``definition.yaml``.

The ``Resistance`` folder is almost exactly the same as the ``Spectra``. Each of the
sources is represented here again (with the same names), and the filename format is the
same, except that the files themselves are all ``.csv``. These measure the resistance
readings of the sources, which are used to derive the physical temperatures of the loads
(against which the spectra are calibrated). Again, each source is allowed to have
multiple "runs" specified by their "run number". However, again in practice it is
very uncommon to have more than one usable run.

The ``S11`` folder contains measurements of the reflection coefficients of the sources,
along with the LNA itself and the internal switch. These are all made with a VNA, and
each reading takes of order a minute. Thus, multiple readings of these measurements can
be taken -- and typically are taken. Inside the ``S11`` folder exist a folder for each
of the main loads (or sources), in which are measurements of the four standards (``open``,
``short``, ``match`` and ``external``). Each of these standards can be measured multiple
times, and so each file has the format ``<standard-name><rep-num>.s1p``, where ``rep-num``
goes from 1 - 99. However, each of the standards for a load is measured one after the
other on the same connection (i.e. there is no disconnection between them, to avoid
issues with different connection characteristics between the standards). Thus, one can't
choose to use repeat number 01 for ``open`` and repeat number 2 for ``short`` for the ``Ambient``
source. For a given source, all standards used must be of the same repeat number (but multiple
runs can *exist* for the source).
Besides the S11's of the sources, we also need measurements of the LNA reflections, and
the internal switch. These exist in the folders ``ReceiverReadingXX`` and ``SwitchingStateXX``
respectively. Here the ``XX`` correspond to what we call a "run" number, which
correspond to a complete re-measuring of the standards at different points in the
observation process. An arbitrary number of these can be performed (up to 99), but only
one is required.

In all cases, the default behaviour of ``edges-io`` is to use the *last* run number and
repeat number available for any given measurement.

Combining Multiple Observations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
As of v0.4.0, ``CalibrationObservation`` objects no longer need to be defined fully by
one directory containing all measurements. While that is still an option (and the easiest
way to define a calibration observation), they can also be defined in a more sophisticated
way internally or externally.

Internally, a ``definition.yaml`` file is allowed (and encouraged) which defines properties of the
observation, and also has ``include`` and ``prefer`` keywords which are used to supplement
or override any particular parts of the observation. For example ``include`` could point
to the top-level of any other observation, which could then be used whenever the
main observation lacks data. If this file exists, by default it is used to construct
the full observation virtually. An incomplete example of such a definition file can
be found `here <example-obs-definition.yaml>`_.


Externally, a different file format is used to explicitly define every single measurement
file in an observation. This is supposed to be exhaustive and complete to make it
unambiguous. An example can be found in the `test-suite <tests/test_data/observation.yaml>`_.
One can use such a file to create a ``CalibrationObservation`` by using the
``CalibrationObservation.from_observation_yaml()`` function.

The way the code actually handles these "virtual" observations is essentially to create
a temporary directory and make symlinks to all the files that are required. This virtual
observation then looks and feels like a normal single observation, but is in fact
patched together from various observations.

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
``HDF5RawSpectrum`` class, which is used to define ``fastspec`` output files.

How the code works in a bit more detail
---------------------------------------
For the sake of developers (lets face it, most users of this particular repo should also
be developers), we will try to explain in a little more detail how the code works here.
This will focus on how the code treats the organization of a calibration observation,
and how it performs checks and makes fixes.

The basic idea is that each directory, and each kind of file, is represented by a
distinct class, describing that kind of thing. For example, the top-level directory
(actually, the top-level plus the ambient temperature directory) is represented by the
``CalibrationObservation`` class, while the ``Spectra`` directory is represented by the
``Spectra`` class, and S1P files are represented by the ``S1P`` class.

All of these classes are subclasses either of ``_DataContainer`` (if it's a folder) or
``_DataFile`` (if it's a file). All of them have a ``path`` attribute which points to its
own path on-disk. ``_DataFile`` classes are much simpler, and typically only know how to
check its own filename for consistency with the specification, and how to read the data
in that particular filetype (they know nothing about their parents).
``_DataContainer`` classes know about their own ``path``,
but also can determine a list of files/subfolders they *contain* (they know nothing about their
parents), and know how to map these files/folders onto their relevant defining classes.
They are able to check their own path for consistency, ensure that all relevant sub-files
exist, ensure that none extra exist, and recursively check the consistency of their
sub-files and folders by calling their checking methods.
Each file and folder in the observation becomes a specific *instance* of one of these
classes (there will be multiple ``S1P`` instances for all of the S11 measurements, and
each may have a different ``name`` attribute to identify the standard it represents).

This top-down hierarchical structure is useful, and similar the to the way Unix filesystems
operate. However, it does mean that a particular instance is not necessarily unique: the
"match" standard S11 will exist within all sources, and since each class doesn't know
its parent, the ``Ambient/Match01.s1p`` cannot be distinguished from the ``HotLoad/Match01.s1p``.
However, a method exists on the top-level ``CalibrationObservation`` which can match a
particular input path to a unique *sequence* of instances which do uniquely define it
(i.e. the first would be a sequence containing a ``LoadS11`` class with ``name=Ambient``
and the second would contain a ``LoadS11`` class with ``name=HotLoad``).

Another thing to note about the setup is the different between the classes and *instances*
of those classes. Much of the functionality of the system is implemented just through
the classes themselves -- one does not need to make instances of the classes to perform
the filesystem checks, for instance. In this case, the ``path`` is given to the ``check()``
method of the class, eg. ``CalibrationObservation.check(path)``, which itself will
call the ``check`` method of any of its children etc. This will never read any data, it
will just check filename formats and contents of directories. However, one can make an
*instance* of the ``CalibrationObservation``, which will itself go and make instances
of all its children, storing them in the top-level class in a nice hierarchical way, in
which each of the children can be used independently. By default, when you create such
an instance, it will first perform the full check that would have been performed (but
in this case it should exit at the first error raised, and raise it as an error, rather
than continuing and printing all errors). Notably, these instances can be used to *read*
the data in the files themselves. The instance will also decide *which* files to use
in the observation (i.e. which run numbers and repeat numbers).



Note
====

This project has been set up using PyScaffold 3.2.3. For details and usage
information on PyScaffold see https://pyscaffold.org/.
