.. intro.rst

.. _maroonxdr_prog_intro:

************
Introduction
************

This Programmer's Manual documents the internal structure of the
MAROON-X DRAGONS pipeline (``maroonxdr``) for developers extending,
maintaining, or porting it. It assumes you have already installed the
development environment and are familiar with the end-to-end reduction
workflow.

For the user-facing instrument description (detector arms, fiber
configurations, frame types, calibration requirements, and known
operational issues), see the *Instrument and Data* chapter of the User Manual.

For installation and the reduction walkthrough, see the Tutorial: the
*Setup and Installation* chapter first, then the CLI or Python API
example chapters for the two versions of the same reduction.

Purpose and Audience
====================

This manual documents the internal structure of ``maroonxdr`` so that a
new developer can extend a primitive, add or modify a recipe, port a
legacy MaroonX algorithm, or track down a bug. It focuses on the parts
of the codebase that DRAGONS itself does not document: the MaroonX tag
system, the recipe-level data flow, the primitive class hierarchy, 
the PDF-report wiring, the parameter ``Config`` classes, and the test 
infrastructure. The goal is to make the design decisions legible so 
that future changes stay consistent with them.

Who this is for:

* Developers adding or modifying primitives, parameters, or recipes.
* Maintainers keeping the pipeline aligned with upstream DRAGONS.
* Anyone porting a legacy MaroonX algorithm into the DRAGONS framework.

.. note::

   MaroonX is currently an *external* DRAGONS instrument package. That
   is why every CLI invocation carries ``--adpkg maroonx_instruments``
   and ``--drpkg maroonxdr``, and why every Python API script begins
   with ``import maroonx_instruments`` followed by
   ``myreduce.drpkg = 'maroonxdr'``. Several design choices covered in
   this manual (the packaging of AstroData tags, the location of
   lookup tables, the calibration database wiring) follow from that
   external-package status. Once upstream integration completes, those
   flags and the explicit import go away, but the internal structure
   documented here will remain.

Architecture Overview
=====================

The MaroonX code is split across two Python packages under ``MAROONXDR/``.
``maroonx_instruments/maroonx/`` defines
the AstroData instrument class ``AstroDataMAROONX``, a set of
``@astro_data_tag``-decorated methods that classify raw files by fiber
configuration, and ``@astro_data_descriptor``-decorated methods that
expose header keywords under stable names. ``maroonxdr/maroonx/`` holds
the primitives, parameter ``Config`` classes, recipes, lookup tables, and
tests, and is where all reduction logic lives. DRAGONS registers the two
packages separately, which is why every CLI invocation passes
``--adpkg maroonx_instruments`` and ``--drpkg maroonxdr``, and why at the 
moment the Python API needs both ``import maroonx_instruments`` and
``myreduce.drpkg = 'maroonxdr'``.

When ``reduce`` (or ``Reduce.runr()``) is called on a file, AstroData opens
it and the ``@astro_data_tag`` methods in ``maroonx_instruments`` assign a
tag set based on the fiber illumination pattern and header keywords. ``recipe_system``
then resolves that tag set to a recipe (the tag-to-recipe mapping is
covered in :ref:`tags`). A recipe is a Python function that instantiates
the ``MAROONX`` primitives class and calls a sequence of primitive methods
on it, passing intermediate products between them.
Each primitive is a method on the ``MAROONX`` class, and its parameter
defaults are drawn from a matching ``Config`` class in a
``parameters_maroonx_*.py`` file. Processed calibrations produced along
the way are registered with ``caldb`` and served back automatically to
later recipes that request them.

The ``MAROONX`` class, defined in ``primitives_maroonx_2D.py``, inherits
from ``Gemini, CCD, NearIR, CalibDBMaroonX``. Those bases contribute the
shared Gemini primitives (``prepare``, ``addDQ``, header standardization,
generic CCD handling) that MaroonX reuses unchanged. On top of them the
class adds the MaroonX-specific primitives: bundle splitting, straylight
removal, stripe extraction, etalon-based wavelength calibration, fiber
combination, and barycentric correction. The full member list is explained
in detail in :ref:`primitives`.

The reduction logic is spread across three primitive files by processing
stage: ``primitives_maroonx_2D.py`` for image-plane operations,
``primitives_maroonx_echelle.py`` for stripe finding and extraction, and
``primitives_maroonx_spectrum.py`` for 1D spectrum-level work including
wavelength assignment and barycentric correction. Recipes live under
``recipes/sq/``, and each ``parameters_maroonx_*.py`` file mirrors the
primitive file of the same suffix.

Repository Layout
=================

The tree below shows the layout of the two MaroonX packages within ``MAROONXDR/``.

.. code-block:: text

   MAROONXDR/
   ├── maroonx_instruments/maroonx/
   │   ├── __init__.py
   │   ├── adclass.py                       # AstroDataMAROONX: tags and descriptors
   │   ├── calibration_maroonx.py           # CalibrationMAROONX: custom caltype resolution
   │   ├── lookup.py                        # header-keyword lookup table
   │   └── tests/
   ├── maroonxdr/maroonx/
   │   ├── primitives_maroonx_2D.py         # image-plane primitives; defines the MAROONX class
   │   ├── primitives_maroonx_echelle.py    # stripe finding, extraction, wavelength calibration
   │   ├── primitives_maroonx_spectrum.py   # 1D spectrum work, fiber combination, barycentric
   │   ├── primitives_calibdb_maroonx.py    # CalibDBMaroonX: custom caldb store/get primitives
   │   ├── parameters_maroonx_2D.py
   │   ├── parameters_maroonx_echelle.py
   │   ├── parameters_maroonx_spectrum.py
   │   ├── parameters_calibdb_maroonx.py
   │   ├── maroonx_utils.py                 # shared helpers used across primitive files
   │   ├── maroonx_plots.py                 # PDF-report plotting helpers
   │   ├── maroonx_echellespectrum/         # subpackage: echelle-spectrum data model
   │   ├── recipes/                         # science-quality and QA recipes
   │   ├── lookups/                         # static calibration tables (BPM, SID, WLS, MDF)
   │   └── tests/
   ├── doc/                                 # this Sphinx documentation tree
   ├── noxfile.py
   ├── pyproject.toml
   └── setup.cfg

Inside ``maroonxdr/maroonx/``, filenames follow a predictable pattern.
Every ``primitives_maroonx_X.py`` file has a sibling
``parameters_maroonx_X.py`` that holds the ``Config`` classes for the
primitives in it; the 2D image primitives, for example, live in
``primitives_maroonx_2D.py`` and their parameters in
``parameters_maroonx_2D.py``.

Recipes live under ``maroonxdr/maroonx/recipes/sq/`` and are named after
the tag set they trigger on: ``recipes_DARK.py`` for darks,
``recipes_BUNDLE.py`` for bundles, and so on. Static calibration data is
grouped by product under ``maroonxdr/maroonx/lookups/`` (``BPM/`` for bad
pixel masks, ``SID/`` for stripe identification, ``WLS/`` for reference
wavelength solutions).

Tests follow the same split as the primitive files: the echelle
primitives are covered by ``maroonxdr/maroonx/tests/echelle_extraction/``,
the 2D primitives by ``tests/image/``, and so on. See :ref:`tests` for
fixtures, markers, and environment variables.

Development Environment
=======================

The full environment setup (cloning the repository, running the ``devenv``
nox session, and configuring the calibration database) is covered in the
*Setup and Installation* chapter of the Tutorial; treat that as the
canonical walkthrough
rather than repeating it here. Day-to-day development happens inside the
``mx_dev`` virtual environment that ``devenv`` creates at ``venv/``, and
that same environment is what the nox sessions below assume is available.
For how the test suite consumes ``DRAGONS_TEST``, see :ref:`tests`.

The developer-facing nox sessions are:

``devenv``
   Creates the ``mx_dev`` virtual environment at ``venv/``, clones or
   updates DRAGONS, and installs ``maroonxdr`` and ``maroonx_instruments``
   in editable mode. Covered in the *Setup and Installation* chapter of the Tutorial; you
   should not need
   to re-run it unless the environment is wiped or dependencies change.

``devconda``
   Conda-based alternative to ``devenv`` that creates the ``mx_devconda``
   environment. Use it only if you specifically need a conda-managed
   stack; the rest of this manual assumes ``mx_dev``.

``preprocess``
   Runs the blessed reduction chain (bundle, dark, flat, wavecal,
   science) that produces the reference data for the regression tests.
   Slow, and only needed by whoever regenerates the blessed data; see
   :ref:`tests`.

``create_inputs``
   Stages the per-module ``inputs/`` and ``refs/`` directories under
   ``DRAGONS_TEST`` from the products of ``preprocess``, by running each
   test module's own staging hook. Also blessing-side only: developers
   who received the packaged test data do not run it.

``unit_tests``
   Runs the synthetic test tier over ``maroonxdr/maroonx/tests/`` and
   ``maroonx_instruments/maroonx/tests/``. Needs no test data. Extra
   pytest arguments are forwarded, so ``nox -s unit_tests -- -k arm``
   will narrow the selection further.

``regression_tests``
   Runs the regression test tier (stored-reference comparison on real
   data, selected with the ``regression`` marker). Needs the test data
   under ``DRAGONS_TEST``.

``coverage``
   Runs both test tiers with ``pytest-cov`` and prints a terminal
   coverage report.

``package_test_data``
   Zips the staged test data into a versioned archive for sharing with
   other developers; see :ref:`tests`.

``package_lookups``
   Zips the lookup FITS files under ``maroonxdr/maroonx/lookups/`` into
   the versioned ``lookups_files_<version>.zip`` that ships with a release;
   see `Building the Lookup Files`_.

``docs``
   Builds this Sphinx documentation using the unified ``doc/conf.py``.
   The ``usermanual``, ``progmanual``, and ``tutorial`` sessions build
   the individual manuals instead, with ``-- --pdf`` for a PDF build.

Run ``nox -l`` from the repository root to see every available session
together with its docstring.

Building the Lookup Files
-------------------------

The lookup FITS files under ``maroonxdr/maroonx/lookups/`` are not tracked
by git. They are built from the legacy HDF5 configuration files and shipped
as ``lookups_files_<version>.zip`` on the GitHub releases page, which users
unzip into ``lookups/`` as described in the Tutorial. Rebuilding them is
only needed when the legacy sources change or when the lookup layout
changes.

Two legacy sources feed the eight runtime files:

.. list-table::
   :header-rows: 1
   :widths: 45 55

   * - Legacy source
     - Lookup files built from it
   * - ``config_b.hdf``, ``config_r.hdf`` (one per arm)
     - ``BPM/BPM_[b|r]_0000.fits``, ``SID/SID_[b|r].fits``,
       ``WLS/WLSTAT_[b|r].fits``
   * - ``wl_combined_final_etalon_peakmodel_2020.hdf`` (shared)
     - ``WLS/REFWAVELENGTH_b.fits``, ``WLS/REFWAVELENGTH_r.fits``

The builders are four functions in ``maroonxdr/maroonx/maroonx_utils.py``.
Each one writes its product through ``astrodata`` and takes the output file
name from the matching lookup dictionary (``maskdb.bpm_dict``,
``siddb.sid_dict``, ``wavelengthdb.statwavelength_dict`` and
``wavelengthdb.refwavelength_dict``), so the registry and the files on disk
cannot drift apart.

``build_bpm_lookup(config_hdf, arm, outdir)``
   Reads ``/bad_pixel_map`` and writes the bad pixel mask as a
   multi-extension file with the mask in a ``SCI`` extension. The legacy
   map uses 1 for good pixels; the lookup follows the DRAGONS DQ
   convention (nonzero is bad), so the builder inverts it.

``build_sid_lookup(config_hdf, arm, outdir)``
   Reads the ``positions`` attribute of ``/identify_stripes`` and writes
   it as the ``SID`` table.

``build_statwls_lookup(config_hdf, arm, outdir)``
   Reads ``/wavelengths_static/fiber_N`` for the five fibers and writes one
   ``FIBER_N`` table per fiber, one column per echelle order.

``build_refwls_lookup(peakmodel_hdf, outdir)``
   Reads the lmfit parameter dump under ``dispersion/parameter`` and the
   ``wls_[blue|red]/fiber_N`` groups and writes both arms in one call: the
   ``PARAMETERS`` table, shared by the two files, and one ``FIBER_N`` table
   per fiber. The legacy parameter dump embeds callables pickled under
   Python 3.7 that current lmfit versions cannot decode, so the builder
   reads the parameter states directly instead of calling
   ``lmfit.Parameters.loads``; no parameter carries an expression, so
   nothing is lost.

``arm`` is ``'BLUE'`` or ``'RED'`` and is written to the ``ARM`` header
card, alongside ``INSTRUME``, ``OBSTYPE``, the legacy file name in
``ORIGNAME`` and the build date in ``DATE``. The functions import ``h5py``
lazily; it is a test-group dependency, so it is available in the ``mx_dev``
environment but is not needed to run the pipeline.

To rebuild the complete set in place, from the repository root in the
``mx_dev`` environment, with the legacy files at their usual location:

.. code-block:: python

   from maroonxdr.maroonx import maroonx_utils as mu

   legacy = '/path/to/legacy/'
   lookups = 'maroonxdr/maroonx/lookups/'
   for arm, letter in (('BLUE', 'b'), ('RED', 'r')):
       config = f'{legacy}config_{letter}.hdf'
       mu.build_bpm_lookup(config, arm, lookups + 'BPM')
       mu.build_sid_lookup(config, arm, lookups + 'SID')
       mu.build_statwls_lookup(config, arm, lookups + 'WLS')
   mu.build_refwls_lookup(
       f'{legacy}wl_combined_final_etalon_peakmodel_2020.hdf',
       lookups + 'WLS',
   )

This overwrites the files under ``lookups/``. Then package the set
for release:

.. code-block:: bash

   nox -s package_lookups

The session zips the eight FITS files and the three READMEs with paths
relative to ``lookups/``, refuses to run if any of them is missing, and
writes ``lookups/lookups_files_<version>.zip`` with the version taken from
``pyproject.toml``.
