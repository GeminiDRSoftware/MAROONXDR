.. _tests:

*****
Tests
*****

Environment Setup
=================

Two environment variables control where test data is found:

``DRAGONS_TEST``
    Root directory for DRAGONS test data (inputs, references,
    preprocessed calibrations). Required for all unit tests::

        export DRAGONS_TEST="/path/to/test/data"

``MAROONX_LEGACY_TEST``
    Root directory containing legacy pipeline reduced data. Required only
    for legacy regression tests::

        export MAROONX_LEGACY_TEST="/path/to/legacy/data"

    The legacy regression tests expect the following structure::

        $MAROONX_LEGACY_TEST/
        └── MaroonX_spectra_reduced/
            ├── Maroonx_masterframes/
            │   └── 202411xx/
            │       ├── darks/               # legacy master darks
            │       └── flats/               # legacy master flats
            └── 20241124/                    # legacy reduced spectra & wavecals


.. note:: The ``devenv`` nox session writes ``DRAGONS_TEST`` into the
   virtualenv activate script automatically. ``MAROONX_LEGACY_TEST`` must
   be set manually if you intend to run legacy regression tests.

Test Data Structure
-------------------

Test data follows the DRAGONS convention. Each test module stores its inputs
under a path derived from the package and test name::

    $DRAGONS_TEST/
    ├── maroonxdr/maroonx/
    │   ├── bundle/
    │   │   ├── test_bundle/inputs/            # raw bundle FITS
    │   │   └── test_bundle_export/inputs/
    │   ├── image/
    │   │   ├── test_file_sorting/inputs/      # debundled FITS
    │   │   ├── test_image_orientation_corrector/inputs/
    │   │   ├── test_ND_filter_check/inputs/
    │   │   ├── test_var/inputs/
    │   │   ├── test_stray_light_removal/inputs/
    │   │   └── test_stripe_finding/inputs/
    │   └── echelle_extraction/
    │       ├── test_extraction/inputs/
    │       ├── test_measure_blaze/inputs/
    │       ├── test_stripe_retrieval/inputs/
    │       └── test_wavecal/inputs/
    ├── maroonx_instruments/maroonx/
    │   ├── test_maroonx/inputs/
    │   └── test_calibration/inputs/
    ├── preprocessed_files/                    # preprocessed calibrations
    │   └── calibrations/
    │       ├── processed_dark/
    │       ├── processed_dark_coeff/
    │       └── processed_flat/
    └── raw_files/                             # raw file download cache

The ``raw_files/`` subdirectory is used as a download cache by the
DRAGONS ``download_from_archive`` helper. The ``inputs/`` directories contain the
actual files each test expects, created by ``create_inputs`` (see below).

The ``preprocessed_files/`` tree holds heavier calibration products (master
darks, flats, dark coefficients) needed by the ``slow`` /
``preprocessed_data`` tests. These are **not** created by the
``create_inputs`` nox session and must be produced separately (e.g. by
running the ``complete_tests`` session or a manual reduction).

Populating Test Data
--------------------

Before running unit tests for the first time, populate the inputs with the
``create_inputs`` nox session::

    nox -s create_inputs

This calls each test module's ``create_inputs()`` function, which downloads
raw files from the Gemini Archive and preprocesses them (e.g. splitting
bundles) into the ``inputs/`` directories listed above.


Running the Tests
=================

Basic Commands
--------------

Run all unit tests::

    pytest maroonxdr/maroonx/tests/ maroonx_instruments/maroonx/tests/

Run a specific test category::

    pytest maroonxdr/maroonx/tests/image/
    pytest maroonxdr/maroonx/tests/echelle_extraction/
    pytest maroonxdr/maroonx/tests/bundle/

Run a single test file::

    pytest maroonxdr/maroonx/tests/image/test_stray_light_removal.py


Test Markers
------------

Five custom markers are registered in ``conftest.py``. Three are
auto-applied based on test location; two are applied manually per test:

.. list-table::
   :header-rows: 1
   :widths: 25 15 60

   * - Marker
     - Applied
     - Description
   * - ``maroonx``
     - auto
     - All tests under ``maroonx/tests/``.
   * - ``regression``
     - auto
     - Tests under any ``regression/`` subdirectory.
   * - ``legacy_regression``
     - auto
     - Tests under ``legacy_regression/``.
   * - ``slow``
     - manual
     - Long-running tests (e.g. full extractions).
   * - ``preprocessed_data``
     - manual
     - Tests that require preprocessed calibrations in ``inputs/``.

Filter tests using markers::

    # Skip slow tests and those needing preprocessed data
    pytest -m "not slow and not preprocessed_data" maroonxdr/maroonx/tests/

    # Run only legacy regression tests
    pytest -m legacy_regression maroonxdr/maroonx/tests/

    # Combine markers
    pytest -m "maroonx and not slow" maroonxdr/maroonx/tests/


Available Tests
===============

The test suite is organized into four main categories, each containing multiple test modules.

**Bundle Processing Tests**

.. autosummary::
   :nosignatures:

   maroonxdr.maroonx.tests.bundle.test_bundle
   maroonxdr.maroonx.tests.bundle.test_bundle_export

**2D Image Processing Tests**

.. autosummary::
   :nosignatures:

   maroonxdr.maroonx.tests.image.test_file_sorting
   maroonxdr.maroonx.tests.image.test_ND_filter_check
   maroonxdr.maroonx.tests.image.test_image_orientation_corrector
   maroonxdr.maroonx.tests.image.test_stray_light_removal
   maroonxdr.maroonx.tests.image.test_stripe_finding
   maroonxdr.maroonx.tests.image.test_var

**Echelle Extraction Tests**

.. autosummary::
   :nosignatures:

   maroonxdr.maroonx.tests.echelle_extraction.test_extraction
   maroonxdr.maroonx.tests.echelle_extraction.test_stripe_retrieval
   maroonxdr.maroonx.tests.echelle_extraction.test_wavecal

**Legacy Regression Tests**

.. autosummary::
   :nosignatures:

   maroonxdr.maroonx.tests.legacy_regression.test_masterdark
   maroonxdr.maroonx.tests.legacy_regression.test_masterflat
   maroonxdr.maroonx.tests.legacy_regression.test_extractions
   maroonxdr.maroonx.tests.legacy_regression.test_fitting
   maroonxdr.maroonx.tests.legacy_regression.test_reduced_science
   maroonxdr.maroonx.tests.legacy_regression.test_reduced_wavecal

GitHub Actions Integration
==========================

The test suite runs automatically on GitHub Actions for every push and pull
request to ``main``, ``release/*``, and ``develop``.

**Workflow file:** ``.github/workflows/testing.yml``

The CI pipeline has two main steps:

1. **Populate test data** — downloads raw files and creates test inputs::

       nox -s create_inputs

2. **Run unit tests** — skipping slow and preprocessed-data tests::

       nox -s unit_tests -- -m "not slow and not preprocessed_data"

Test data is cached between runs (``~/mx_test``) so that archive downloads
are only performed once.


Missing or Desirable Tests
==========================

.. todo::

   Document missing or desirable tests for MAROONX data reduction pipeline.

   maroonx/maroonx_fit/
   maroonx/maroonx_echellespectrum/

