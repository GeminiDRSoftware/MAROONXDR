.. _tests:

*****
Tests
*****

Environment Setup
=================

Two environment variables control where the test suite finds its data.

``DRAGONS_TEST``
   Root of the DRAGONS test tree (per-test ``inputs/``, raw file cache,
   preprocessed calibrations). Required for every test session. The
   ``devenv`` nox session appends an ``export`` for this variable to
   ``venv/bin/activate``, so simply activating ``mx_dev`` (see
   :ref:`maroonx_setup`) sets it for you. To override the location, or to
   set it outside a ``devenv``-managed shell, export it manually:

   .. code-block:: bash

      export DRAGONS_TEST="/path/to/mx_test"

``MAROONX_LEGACY_TEST``
   Root of the legacy MaroonX data used by the ``legacy_regression``
   suite. Contains two top-level subdirectories: ``MaroonX_spectra/``
   (raw legacy frames) and ``MaroonX_spectra_reduced/`` (legacy pipeline
   configuration files, master frames, and wavelength solutions).
   Optional, and only consulted by tests under
   ``maroonxdr/maroonx/tests/legacy_regression/``. Not set by ``devenv``;
   export it manually if you intend to run those tests:

   .. code-block:: bash

      export MAROONX_LEGACY_TEST="/path/to/legacy/data"

The sections that follow assume both variables are set
correctly. For the full list of nox sessions, see the Programmer Manual 
introduction: :ref:`maroonxdr_prog_intro`.

Test Data Layout
================

DRAGONS test data follows a convention where each test module's inputs
live under a path derived from that module's location. The layout below
shows the top level of ``$DRAGONS_TEST`` and where each test group's
data lives inside it.

.. code-block:: text

   $DRAGONS_TEST/
   ├── raw_files/                          # raw FITS bundle download cache
   ├── preprocessed_files/
   │   └── calibrations/
   │       ├── processed_dark/             # master darks per exptime x arm
   │       ├── processed_dark_coeff/       # dark scaling coefficients per arm
   │       ├── processed_flat/             # master flats per arm
   │       └── processed_wavecal/          # dynamic etalon wavelength solutions
   ├── maroonxdr/maroonx/
   │   ├── bundle/<test_name>/inputs/            # bundle-splitting tests
   │   ├── image/<test_name>/inputs/             # 2D image-plane tests
   │   ├── echelle_extraction/<test_name>/inputs/  # stripe and extraction tests
   │   └── legacy_regression/<test_name>/inputs/ # legacy pipeline comparison
   └── maroonx_instruments/maroonx/
       ├── test_maroonx/inputs/
       └── test_calibration/inputs/

The pytest-dragons plugin resolves each test module's data by mapping
the module's on-disk path to a matching path under ``$DRAGONS_TEST``.
For example, the test module at
``maroonxdr/maroonx/tests/image/test_stripe_finding.py`` finds its
inputs under ``$DRAGONS_TEST/maroonxdr/maroonx/image/test_stripe_finding/inputs/``.
Regression tests that also carry stored reference outputs follow the
same convention with a sibling ``refs/`` directory.

The ``raw_files/`` subdirectory is a shared download cache populated by
the ``download_raws`` nox session, which pulls raw bundles from the
Gemini Archive. Per-module ``create_inputs()`` hooks read from this
cache when preparing their ``inputs/`` directories, so it does not
follow the per-module convention itself.

The ``preprocessed_files/`` tree is a shared calibration store
populated by the ``preprocess`` nox session, which runs full
reductions (bundle, dark, flat, wavecal, science) and lands the
processed calibrations under ``calibrations/<caltype>/``. Tests marked
``@pytest.mark.preprocessed_data`` read from this tree. Like
``raw_files/``, it sits outside the per-module layout.

Populating Test Data
====================

Three nox sessions populate test data. Run ``download_raws`` once to
bootstrap the raw file cache, then use ``create_inputs`` or ``preprocess``
depending on which tests you plan to run.

``download_raws``
   Fills the shared raw file cache at ``$DRAGONS_TEST/raw_files/`` by
   downloading every file listed in the ``MAROONX_TEST_MANIFEST`` dict in
   ``maroonxdr/maroonx/tests/conftest.py`` from the Gemini Observatory
   Archive. Files already present are skipped. If the archive returns
   HTTP 403 for a file (proprietary data), the session emits a warning
   and continues, so downstream tests may fail later with missing files.
   Rerun after the manifest changes or on a fresh ``DRAGONS_TEST`` root.

``create_inputs``
   Runs each test module's ``create_inputs()`` hook, which reads from
   ``raw_files/`` and writes per-module ``inputs/`` directories under
   ``$DRAGONS_TEST``. Required before the first ``unit_tests`` run. Not
   every test module has a hook; those that don't have it rely on fixtures 
   or the shared calibration store.

``preprocess``
   Runs the five preprocessing scripts under
   ``maroonxdr/maroonx/tests/preprocess/`` in order (``bundle``, ``dark``,
   ``flat``, ``wavecal``, ``science``). This step is slow. Outputs land in
   the shared calibration store at ``$DRAGONS_TEST/preprocessed_files/calibrations/``
   and in the ``legacy_regression/`` per-module ``inputs/`` directories.


Pytest Markers
==============

Five custom markers are registered in ``conftest.py``. Three are attached
automatically by ``pytest_collection_modifyitems`` based on a substring
match against the test file's path, and two are applied by hand with a
``@pytest.mark`` decorator on the test.

``maroonx`` (auto)
   Attached to every collected item whose path contains
   ``maroonx/tests``, which covers the entire suite. Use it to scope a
   run when other test roots are on the command line, for example
   ``pytest -m maroonx maroonxdr/ maroonx_instruments/``.

``regression`` (auto)
   Attached to every item whose path contains the substring
   ``regression``. Because ``legacy_regression`` also matches, legacy
   tests carry both this marker and ``legacy_regression``. Example:
   ``pytest -m regression maroonxdr/maroonx/tests/``.

``legacy_regression`` (auto)
   Attached to every item whose path contains ``legacy_regression``,
   which are the tests that compare DRAGONS outputs against frozen
   legacy pipeline products. Requires ``MAROONX_LEGACY_TEST``. Example:
   ``pytest -m legacy_regression maroonxdr/maroonx/tests/``.
   
   .. note::

      Legacy regression is not directly achievable because of ``photutils.Background2D`` changes in v2.0.0.
      The legacy pipeline used v0.7.1, which is no longer available in the current environment.
      These tests justify the existence of the ``removeStrayLight_legacyPatch`` primitive that
      feeds from custom intermediate legacy outputs, ``.npy`` files.

``slow`` (manual)
   Decorate long-running tests (for example full extractions) with
   ``@pytest.mark.slow``. The CI workflow deselects them with
   ``pytest -m "not slow" maroonxdr/maroonx/tests/``.

``preprocessed_data`` (manual)
   Decorate tests that read from the shared
   ``$DRAGONS_TEST/preprocessed_files/`` calibration store with
   ``@pytest.mark.preprocessed_data``. Skip them with
   ``pytest -m "not preprocessed_data" maroonxdr/maroonx/tests/``.


Fixtures and Test Utilities
===========================

``conftest.py`` provides the following fixtures and one utility for use
in new tests.

Session-scoped fixtures
-----------------------

``dragons_test_root``
   Returns the ``$DRAGONS_TEST`` root as a ``pathlib.Path``. Skips the
   test if the environment variable is unset or the path does not exist.

``legacy_test_root``
   Returns the ``$MAROONX_LEGACY_TEST`` root as a ``pathlib.Path``. Skips
   the test if the environment variable is unset or the path does not
   exist.

``preprocessed_files_path``
   Returns ``$DRAGONS_TEST/preprocessed_files``. Skips if the directory
   is missing.

Function-scoped path fixtures
-----------------------------

``processed_dark_path``
   Returns ``$DRAGONS_TEST/preprocessed_files/calibrations/processed_dark``.
   Skips if the directory is missing.

``processed_dark_coeff_path``
   Returns ``$DRAGONS_TEST/preprocessed_files/calibrations/processed_dark_coeff``.
   Skips if the directory is missing.

``path_to_legacy_darks``, ``path_to_legacy_flats``, ``path_to_legacy_wavecal``, ``path_to_legacy_science``, ``path_to_legacy_reduced``, ``path_to_legacy_bkg``
   Return per-caltype subdirectories under ``$MAROONX_LEGACY_TEST``
   holding the legacy pipeline products used by
   ``legacy_regression`` tests. Each skips if its resolved directory is
   missing. ``path_to_legacy_bkg`` is session-scoped; the other five are
   function-scoped.

Synthetic AstroData fixtures
----------------------------

``ad_min``
   Yields a minimal MaroonX ``AstroData`` object backed by a 4400x4400
   ones array, with the headers required for tag resolution and
   descriptor lookup. Parameterized over ``'RED'`` and ``'BLUE'``.

``ad_echelle``
   Yields an ``AstroData`` object with synthetic ``STRIPES``,
   ``F_STRIPES``, and ``STRIPES_MASKS`` sparse-matrix dicts built from
   the SID lookup so echelle-stage primitives can run without real
   extracted data. Parameterized over ``'RED'`` and ``'BLUE'``.

``arm``
   Yields the string ``'BLUE'`` or ``'RED'``, for tests that only need
   the arm label rather than a full ``AstroData`` object.

Utility
-------

``assert_allclose_with_max_fails``
   Compares two arrays with ``rtol``/``atol`` tolerances and a
   ``max_fails`` budget: raises ``AssertionError`` only if more than
   ``max_fails`` elements are outside tolerance, and (by default)
   triggers ``pytest.xfail`` when the number of failing elements is
   nonzero but within budget. Use this in regression tests where a small
   number of pixel-level differences are acceptable.


Running the Tests
=================

For the general nox-sessions review see :ref:`maroonxdr_prog_intro`; this
section covers what each test-running session does and how to pass
pytest arguments through.

``unit_tests``
   Runs the fast pytest suite over ``maroonxdr/maroonx/tests/`` and
   ``maroonx_instruments/maroonx/tests/``, excluding
   ``legacy_regression/``. Assumes ``create_inputs`` has already run.

``regression_tests``
   Currently runs the ``regression/`` subdirectory under
   ``maroonxdr/maroonx/tests/``. This subdirectory contains only
   ``__init__.py`` right now (no tests), so the session runs but
   collects nothing until regression tests are added.

``legacy_regression_tests``
   Runs the comparisons under
   ``maroonxdr/maroonx/tests/legacy_regression/``. Assumes ``preprocess
   -- --legacy-patch`` has already run and that ``MAROONX_LEGACY_TEST``
   is set.

``coverage``
   Same scope as ``unit_tests`` but with ``--cov`` reporting.

Any argument after ``--`` in ``nox -s <session> -- <args>`` is forwarded
to pytest. For example, to skip slow and preprocessed-data tests in a
unit-test run:

.. code-block:: bash

   nox -s unit_tests -- -m "not slow and not preprocessed_data"

This is the same passthrough pattern that CI uses.


Available Tests
===============

The test suite is organized under ``maroonxdr/maroonx/tests/`` (pipeline tests)
and ``maroonx_instruments/maroonx/tests/`` (instrument-package tests), each
grouped by processing stage or subject.

Pipeline tests (``maroonxdr/maroonx/tests/``)
---------------------------------------------

**bundle/** (bundle splitting and reduced-bundle export)

* ``test_bundle.py``
* ``test_bundle_export.py``

**image/** (2D image processing)

* ``test_file_sorting.py``
* ``test_image_orientation_corrector.py``
* ``test_ND_filter_check.py``
* ``test_stack.py``
* ``test_stray_light_removal.py``
* ``test_stripe_finding.py``
* ``test_subtract_overscan.py``
* ``test_var.py``

**echelle_extraction/** (stripe retrieval, extraction, wavelength calibration)

* ``test_extraction.py``
* ``test_measure_blaze.py``
* ``test_stripe_retrieval.py``
* ``test_synthetic_dark.py``
* ``test_wavecal.py``

**legacy_regression/** (comparison against the legacy pipeline; the
``legacy_adapter.py`` module in this directory is a helper, not a test module)

* ``test_extractions.py``
* ``test_fitting.py``
* ``test_masterdark.py``
* ``test_masterflat.py``
* ``test_reduced_science.py``
* ``test_reduced_wavecal.py``

.. note::

   The ``regression/`` subdirectory exists as an empty stub (only
   ``__init__.py``) and is not yet populated.

Instrument-package tests (``maroonx_instruments/maroonx/tests/``)
-----------------------------------------------------------------

Tests that exercise the AstroData instrument definition and calibration
descriptors.

* ``test_calibration.py``
* ``test_maroonx.py``

GitHub Actions Integration
==========================

The workflow file ``.github/workflows/testing.yml`` runs on every push and
pull request to ``main``, ``release/*``, and ``develop``, and can also be
triggered manually via ``workflow_dispatch``.

Pipeline steps, in order:

* Restore (or seed) the ``~/mx_test`` cache. The cache key includes a
  ``CACHE_NUMBER`` env var; bumping it invalidates the cache and forces
  a fresh download on the next run.
* Set up Python 3.12.
* Export ``DRAGONS_TEST=$HOME/mx_test`` into the job environment.
* Install ``nox``.
* Run ``nox -s download_raws-3.12`` to pull raw files from the archive.
* Run ``nox -s create_inputs-3.12`` to build per-test input directories.
* Run ``nox -s unit_tests-3.12 -- -m "not slow and not preprocessed_data"``.

What CI does *not* run: the ``preprocess`` and ``legacy_regression_tests``
sessions. In practice this means any test marked ``preprocessed_data`` and
every test under ``legacy_regression/`` is skipped in CI, so a green CI run
does not prove that the legacy comparisons still pass. Developers should
run those sessions locally before merging changes to the reduction path.


Missing or Desirable Tests
==========================

The following areas have no unit test coverage yet:

* ``maroonx/maroonx_fit/``
* ``maroonx/maroonx_echellespectrum/``


