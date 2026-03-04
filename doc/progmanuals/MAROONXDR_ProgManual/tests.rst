.. _tests:

*****
Tests
*****


Set the required environment variable::

    export DRAGONS_TEST="/path/to/test/data"

In this path will be created a ``science_dir/`` directory and automatically download test files from the Gemini Archive on first run if they are not already present.


Running the Tests
=================

Basic Commands
--------------

Run all tests::

    pytest maroonxdr/maroonx/tests/

Run specific test category::

    pytest maroonxdr/maroonx/tests/image/
    pytest maroonxdr/maroonx/tests/regression/
    pytest maroonxdr/maroonx/tests/echelle_extraction/
    pytest maroonxdr/maroonx/tests/bundle/

Run a single test file::

    pytest maroonxdr/maroonx/tests/image/test_stray_light_removal.py

Verbose output for debugging::

    pytest -v -s maroonxdr/maroonx/tests/

Test Markers
------------

The test suite defines three custom markers that are auto-applied based on test location:

- ``maroonx``: Applied to all tests in ``maroonx/tests/``
- ``regression``: Applied to all tests in ``regression/`` subdirectories
- ``slow``: Must be manually applied to slow-running tests

Filter tests using markers::

    # Skip slow tests
    pytest -m "not slow" maroonxdr/maroonx/tests/

    # Skip regression tests
    pytest -m "not regression" maroonxdr/maroonx/tests/

    # Run only regression tests
    pytest -m regression maroonxdr/maroonx/tests/

    # Combine markers
    pytest -m "maroonx and not slow" maroonxdr/maroonx/tests/

Command-Line Flags
------------------

``--preprocess-bundles``
    Runs bundle preprocessing before tests. This splits bundle files containing both
    arms into separate BLUE and RED FITS files. Required for tests that depend on
    split bundle data::

        pytest --preprocess-bundles maroonxdr/maroonx/tests/


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

**Regression Tests**

.. autosummary::
   :nosignatures:

   maroonxdr.maroonx.tests.regression.test_masterdark
   maroonxdr.maroonx.tests.regression.test_masterflat
   maroonxdr.maroonx.tests.regression.test_extractions
   maroonxdr.maroonx.tests.regression.test_fitting
   maroonxdr.maroonx.tests.regression.test_reduced_science

GitHub Actions Integration
==========================

The test suite runs automatically on GitHub Actions for every push and pull request.

**Workflow file:** ``.github/workflows/testing.yml``

**Running tests:**

Tests are executed via nox, skipping slow tests and preprocessing bundles::

    nox -s unit_tests-3.12 -- -m "not slow" --preprocess-bundles

This runs:
  - The ``unit_tests`` nox session
  - Only non-slow tests (``-m "not slow"``)
  - Automatically download the test manifest and splits the bundles (``--preprocess-bundles``)


Missing or Desirable Tests
==========================

.. todo::

   Document missing or desirable tests for MAROONX data reduction pipeline.

   maroonx/maroonx_fit/
   maroonx/maroonx_echellespectrum/

