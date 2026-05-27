.. maroonx_reduction_api.rst

.. _maroonx_reduction:

****************************************
MAROON-X DRP: Using the Reduce class API
****************************************

Complete Reduction Workflow Using the Python API
=================================================

This guide walks through the end-to-end reduction of a single MAROON-X
science exposure using the DRAGONS Python API (the ``Reduce`` class and
``dataselect.select_data``). The example reduces an exposure of object
HD 3651 taken on 2025-07-17. The final output is a fully reduced,
wavelength-calibrated spectrum with science quality.

The same workflow expressed through the DRAGONS command-line tools is in
:ref:`maroonx_reduction_cli`. Both tutorials reduce the same files and
produce the same outputs; they differ only in interface.

Overview of Reduction Steps
---------------------------

The full reduction proceeds in nine steps, executed in order:

1. **Debundle** - split raw GOA bundles into per-arm FITS files
2. **Master Darks** - combine raw darks per exposure time and arm
3. **Dark Coefficients** - fit dark scaling vs. exposure time per arm
4. **Master Flats** - build per-arm master flats from DFFFD + DDDDF inputs
5. **Wavelength Calibration** - fit the dynamic etalon wavelength solution
6. **Synthetic Darks** - interpolate a dark matched to each science exposure
7. **Science Reduction** - extract and wavelength-calibrate science spectra
8. **Barycentric Correction** - compute the BERV shift of the reduced spectra
9. **Export Bundle** - combine BLUE and RED arms into a single output FITS

.. note:: MAROON-X has two independent arms (BLUE and RED). Every
   calibration and science step is run separately per arm.

.. important:: MAROON-X is currently an external instrument package, so
   every reduction script must (a) ``import maroonx_instruments`` to
   register the AstroData tags, and (b) set ``myreduce.drpkg = 'maroonxdr'``
   on every ``Reduce`` instance to load the maroonxdr recipe library.
   These two lines will no longer be needed once MAROON-X is integrated
   into DRAGONS.

Prerequisites
=============

Before beginning the reduction, ensure you have:

* A working ``maroonxdr`` installation in a Python virtual environment, with
  the DRAGONS framework available on the same environment. See the
  installation guide of the User Manual for details.
* All raw FITS files for the example (darks, flats, wavecal, science)
  downloaded from the Gemini Observatory Archive (GOA) into a single
  working directory. We call this directory ``science_dir/`` throughout the
  tutorial.

Directory Structure
-------------------

The tutorial assumes a single working directory holding every raw bundle::

    science_dir/
    └── N*.fits                   # Raw GOA bundle files

As the reduction proceeds, ``Reduce.runr()`` writes its outputs into the
working directory and registers calibrations under a ``calibrations/``
subtree managed by ``caldb``::

    science_dir/
    ├── N*.fits                          # Raw GOA bundle files
    ├── *_b_*.fits, *_r_*.fits           # Debundled per-arm files and
    │                                    # downstream reduction products
    ├── calibrations/
    │   ├── processed_dark/              # Master darks (per exptime, per arm)
    │   ├── processed_dark_coeff/        # Dark scaling coefficients (per arm)
    │   ├── processed_flat/              # Master flats (per arm)
    │   └── processed_wavecal/           # Dynamic etalon wavelength solutions
    └── reduce_*.log                     # Reduction log files

.. note:: Calibration files produced by ``Reduce.runr()`` are written
   twice: once in ``science_dir/`` and once in the matching
   ``calibrations/<caltype>/`` subdirectory. The copy under
   ``calibrations/`` is what ``caldb`` indexes and serves to later steps.

.. note:: This tutorial assumes ``caldb`` is configured and initialised.
   See :ref:`maroonx_caldb_setup` for a one-time setup walkthrough.

Dataset
-------

This example reduces a single 300s exposure of HD3651 taken on 2025-07-17, 
against calibrations taken on 2025-07-01, 2025-07-17, and 2025-07-21.

All files can be downloaded from the Gemini Observatory Archive (GOA). The
table below lists the raw GOA bundle ranges; each ``N…M*.fits`` bundle
contains one BLUE and one RED exposure.

+----------+-------------------------------+----------------------------------------------+
| Type     | GOA bundle range              | Notes                                        |
+==========+===============================+==============================================+
| Darks    | ``N20250717M6066`` -          | DDDDE pattern, BLUE and RED, at exposure     |
|          | ``N20250721M6059``            | times 60, 120, 300, 600, 900, 1200, 1800 s   |
+----------+-------------------------------+----------------------------------------------+
| Flats    | ``N20250701M6126`` -          | DFFFD and DDDDF illumination patterns,       |
|          | ``N20250701M6271``            | BLUE and RED                                 |
+----------+-------------------------------+----------------------------------------------+
| Wavecal  | ``N20250717M5948``            | DEEEE etalon frame, BLUE and RED             |
+----------+-------------------------------+----------------------------------------------+
| Science  | ``N20250717M5299``            | HD 3651, 300 s, SOOOE pattern,               |
|          |                               | BLUE and RED                                 |
+----------+-------------------------------+----------------------------------------------+

The dark exposure times cover the range needed by the dark-coefficient fit
(Step 3), which interpolates a synthetic dark matched to the science
exposure time.




Reduction Steps
===============


Step 0: Environment Setup
--------------------------

Import the libraries used throughout the tutorial, configure logging, and
define a small ``get_files`` helper plus the per-arm and per-exposure-time
tag lists that the loops in later steps iterate over.

.. code-block:: python

    import os
    import itertools as it
    from pathlib import Path

    from gempy.adlibrary import dataselect
    from gempy.utils import logutils
    from recipe_system.reduction.coreReduce import Reduce

    import astrodata
    import maroonx_instruments  # noqa - registers MAROON-X AstroData tags

    # Configure logging. Use mode='debug' for verbose output when troubleshooting.
    logutils.config(file_name='reduce.log', mode='standard')

    # Working directory. Change the path here to suit your installation.
    science_dir = Path('/path/to/science_dir')
    os.chdir(science_dir)

    def get_files(pattern='*.fits'):
        """Return a sorted list of files in science_dir matching pattern."""
        return sorted(str(p) for p in science_dir.glob(pattern))

    # Tag lists used by the per-arm and per-exposure-time loops below
    exptime_tags = ['60s', '120s', '300s', '600s', '900s', '1200s', '1800s']
    arm_tags = ['BLUE', 'RED']

.. note:: ``get_files`` is the Python analogue of the CLI's ``*.fits`` glob.
   Later steps reuse it with narrower patterns (``*_reduced.fits``,
   ``*_barycor.fits``) to mirror the CLI's targeted ``dataselect`` calls.



Step 1: Debundle Raw Data
-------------------------

**Purpose**: split each raw GOA bundle into one BLUE and one RED FITS file.

MAROON-X data arrives from GOA as bundle files (tagged ``BUNDLE``) that
contain both arms in a single FITS. The default ``processBundle`` recipe
picks up every ``BUNDLE``-tagged file and writes a per-arm FITS for each
one. For example, the science bundle ``N20250717M5299.fits`` yields
``20250717T144308Z_SOOOE_b_0300.fits`` and
``20250717T144308Z_SOOOE_r_0004.fits``.

Select every bundle in the working directory and run ``Reduce``:

.. code-block:: python

    # Select bundles
    selected_bundles = dataselect.select_data(get_files(), tags=['BUNDLE'])

    # Debundle (no recipename set: processBundle is the default for BUNDLE)
    myreduce = Reduce()
    myreduce.files.extend(selected_bundles)
    myreduce.drpkg = 'maroonxdr'
    myreduce.runr()

After this step ``science_dir/`` contains the original ``N*.fits`` bundles
plus the per-arm debundled files used by every subsequent step.


Step 2: Master Darks
---------------------

**Purpose**: combine the raw darks at each exposure time into a master
dark, per arm.

For every (exposure time, arm) pair, ``dataselect.select_data`` finds the
matching raw darks (DDDDE pattern, tagged ``RAW,DARK``) and ``Reduce``
runs the default ``makeProcessedDark`` recipe. The dataset has 7 exposure
times × 2 arms, so 14 master darks total will be created.

For a single combination of exposure time and arm, the process is:

.. code-block:: python

    # Process darks for a specific exposure time and arm, 300s and BLUE
    selected_darks = dataselect.select_data(get_files(), tags=['RAW', 'DARK', '300s', 'BLUE'])

    myreduce = Reduce()
    myreduce.files.extend(selected_darks)
    myreduce.drpkg = 'maroonxdr'
    myreduce.runr()

Run the full loop over every exposure time and arm:

.. code-block:: python

    # Process darks for each exposure time and arm
    for exptime, arm in it.product(exptime_tags, arm_tags):

        selected_darks = dataselect.select_data(
            get_files(), tags=['RAW', 'DARK', exptime, arm])

        myreduce = Reduce()
        myreduce.files.extend(selected_darks)
        myreduce.drpkg = 'maroonxdr'
        myreduce.runr()

Each output is a single master dark file named after the first raw frame
in the stack with a ``_dark`` suffix - for example
``20250721T170049Z_DDDDE_b_0300_dark.fits``. The master darks land in
``calibrations/processed_dark/`` (and as duplicates in ``science_dir/``).

**Verify processed darks:**

.. code-block:: python

    # List all processed darks
    print(dataselect.select_data(get_files(), tags=['PROCESSED', 'DARK']))

Note the presence of the ``PROCESSED`` tag, which serves to distinguish the processed
darks from the raw darks that are present in the directory.

.. note::
    Processed calibration files are saved twice: in the current working
    directory ``science_dir/`` and in the matching ``calibrations/<caltype>/``
    subdirectory. The copy under ``calibrations/`` is the one ``caldb``
    indexes and serves to later steps.

Step 3: Dark Coefficients
--------------------------

**Purpose**: fit, per arm, the dark scaling as a function of exposure time
so a synthetic dark can later be interpolated for any science exposure.

The ``makeDarkCoefficients`` recipe takes every master dark for one arm
(one per exposure time) and produces a single dark-coefficient file. The
input selector takes all ``PROCESSED,DARK`` files for the arm from the
working directory and excludes the ``DARK_COEFF`` and ``DARK_SYNTH`` tags
so previous coefficient or synthetic-dark outputs in the same directory
are not re-ingested.

Run once per arm:

.. code-block:: python

    for arm in arm_tags:

        # Select all master darks for this arm, excluding any previous
        # coefficient or synthetic-dark outputs
        selected_darks = dataselect.select_data(
            get_files(),
            tags=['PROCESSED', 'DARK', arm],
            xtags=['DARK_COEFF', 'DARK_SYNTH'])

        # Fit the dark-coefficient model
        myreduce = Reduce()
        myreduce.files.extend(selected_darks)
        myreduce.drpkg = 'maroonxdr'
        myreduce.recipename = 'makeDarkCoefficients'
        myreduce.runr()

Each call writes one ``*_darkCoefficients.fits`` file (named after the
first input master dark) into ``science_dir/``, and a copy under the
``processed_dark_coeff`` caltype in ``calibrations/processed_dark_coeff/``.
The output carries the ``DARK_COEFF`` tag, which is how later steps locate
it.

**Verify dark-coefficient files:**

.. code-block:: python

    # List dark coefficient files
    print(dataselect.select_data(get_files(), tags=['DARK_COEFF']))


Step 4: Master Flats
---------------------

**Purpose**: build a master flat per arm.

MAROON-X flats come in two illumination patterns: ``DFFFD`` (flat on
fibers 2-4) and ``DDDDF`` (flat on fiber 5 only). The dataset for this
tutorial contains both patterns and uses the ``makeProcessedFlatDFFFF``
recipe, which combines the two into a single ``DFFFF`` master flat per
arm.

Run the reduction per arm:

.. code-block:: python

    for arm in arm_tags:

        # Select every raw flat for this arm (both DFFFD and DDDDF)
        selected_flats = dataselect.select_data(get_files(), tags=['RAW', 'FLAT', arm])

        # Build the master flat
        myreduce = Reduce()
        myreduce.files.extend(selected_flats)
        myreduce.drpkg = 'maroonxdr'
        myreduce.recipename = 'makeProcessedFlatDFFFF'
        myreduce.runr()

Each call writes one ``*_flat.fits`` master flat into ``science_dir/`` and
a copy under the ``processed_flat`` caltype in
``calibrations/processed_flat/``. The output carries the ``PROCESSED,FLAT``
tag set.

.. note:: The default recipe ``makeProcessedFlat`` combines ``DFFFD`` and
   ``FDDDF`` patterns into an ``FFFFF`` master flat. When the dataset has
   ``DDDDF`` instead of ``FDDDF`` (as in this tutorial), use the
   ``makeProcessedFlatDFFFF`` variant explicitly via ``recipename``.

**Verify processed flats:**

.. code-block:: python

    # List all processed flats
    print(dataselect.select_data(get_files(), tags=['PROCESSED', 'FLAT']))


Step 5: Wavelength Calibration
-------------------------------

**Purpose**: fit, per arm, a dynamic wavelength solution from the etalon
calibration frame.

The Fabry-Perot etalon is the primary wavelength calibration file for MAROON-X.
Anchored to an initial ThAr solution, the etalon provides the refined
wavelength and drift solution that the science reduction (Step 7) applies
to extracted spectra.

The ``makeDynamicWavecal`` recipe takes the raw etalon frame
(``DEEEE`` pattern, tagged ``RAW,ETALON``) for one arm, identifies the
peaks in every order, and fits the per-pixel wavelength solution. The
``getPeaksAndPolynomials:multithreading=True`` override parallelizes the
fit and is recommended on multi-core machines.

Run once per arm:

.. code-block:: python

    for arm in arm_tags:

        # Select the raw etalon wavecal frame for this arm
        selected_wavecal = dataselect.select_data(get_files(), tags=['RAW', 'ETALON', arm])

        # Fit the dynamic wavelength solution
        myreduce = Reduce()
        myreduce.files.extend(selected_wavecal)
        myreduce.drpkg = 'maroonxdr'
        myreduce.recipename = 'makeDynamicWavecal'
        myreduce.uparms = {'getPeaksAndPolynomials:multithreading': True}
        myreduce.runr()

Each call writes one ``*_wavecal.fits`` file into ``science_dir/`` and a
copy under the ``processed_wavecal`` caltype in
``calibrations/processed_wavecal/``. The output carries the
``PROCESSED,WAVECAL`` tag set.

**Verify wavelength calibrations:**

.. code-block:: python

    # List processed wavecal files
    print(dataselect.select_data(get_files(), tags=['PROCESSED', 'WAVECAL']))


Step 6: Synthetic Darks
------------------------

**Purpose**: interpolate, per science frame, a synthetic dark matched to
its exact exposure time using the dark coefficients fit in Step 3.

The ``makeSyntheticDark`` recipe takes the raw science frames and uses the
``DARK_COEFF`` calibration to build a dark frame at the science exposure
time. The synthetic darks are stored under the ``processed_dark`` caltype
with the ``DARK_SYNTH`` tag.

Run once per arm:

.. code-block:: python

    for arm in arm_tags:

        # Select raw science frames for this arm
        selected_sci = dataselect.select_data(get_files(), tags=['RAW', 'SCI', arm])

        # Build the synthetic dark for each science frame
        myreduce = Reduce()
        myreduce.files.extend(selected_sci)
        myreduce.drpkg = 'maroonxdr'
        myreduce.recipename = 'makeSyntheticDark'
        myreduce.runr()


Step 7: Science Data Reduction
-------------------------------

**Purpose**: extract and calibrate the science spectra, applying
the master flat, the dynamic wavelength solution, and the synthetic dark
matched to the science exposure time.

The default recipe for ``RAW,SCI`` files runs the full echelle reduction:
dark and flat correction, straylight removal, stripe extraction, optimal
extraction, wavelength assignment, and fiber combination. ``caldb``
resolves each input calibration automatically from the work done in
Steps 2-6.

The three ``uparms`` overrides below are the recommended defaults for the
current pipeline: ``extractStripes:straylight_removal_fibers=[5]`` applies
the straylight correction to the calibration fiber, ``getPeaksAndPolynomials``
runs the per-order fits in parallel, and ``combineFibers:max_clips=20``
sets the outlier rejection threshold when combining the science fibers.

Run once per arm:

.. code-block:: python

    for arm in arm_tags:

        # Select raw science frames for this arm
        selected_sci = dataselect.select_data(get_files(), tags=['RAW', 'SCI', arm])

        # Full echelle reduction
        myreduce = Reduce()
        myreduce.files.extend(selected_sci)
        myreduce.drpkg = 'maroonxdr'
        myreduce.uparms = {
            'extractStripes:straylight_removal_fibers': [5],
            'getPeaksAndPolynomials:multithreading': True,
            'combineFibers:max_clips': 20,
        }
        myreduce.runr()

Each call writes a ``*_reduced.fits`` file per input science frame into
``science_dir/``. The reduced products carry the ``PROCESSED_SCIENCE`` tag
and are the input to the barycentric-correction step (Step 8).

**Verify reduced science frames:**

.. code-block:: python

    # List all reduced science frames
    print(dataselect.select_data(get_files(), tags=['PROCESSED_SCIENCE']))


Step 8: Barycentric Correction
-------------------------------

**Purpose**: compute the solar-system barycentric correction.

The ``applyBarycentricCorrection`` recipe queries SIMBAD to
resolve the target's coordinates, computes the barycentric Earth radial
velocity (BERV), and writes it to ``*_barycor.fits``.

The ``barycentricCorrection`` primitive accepts two related parameters:

* ``target_name`` — substring match against the FITS ``OBJECT`` header. When
  set, only files whose ``OBJECT`` contains this string are processed; when
  ``None``, every input is processed using its own ``OBJECT`` value for
  the SIMBAD lookup.
* ``simbad_target_name`` — explicit SIMBAD-resolvable name. Use when the
  FITS ``OBJECT`` value (e.g. ``HD3651``) is not directly resolvable and
  you need to override the lookup string (e.g. ``HD 3651``).

For this tutorial the FITS header has ``OBJECT='HD3651'`` (no space), which
is already SIMBAD-resolvable, so only ``target_name`` is needed.

Run once per arm:

.. code-block:: python

    for arm in arm_tags:

        # Select the reduced science file for this arm
        selected_reduced = dataselect.select_data(
            get_files('*_reduced.fits'), tags=['PROCESSED_SCIENCE', arm])

        # Apply the barycentric correction
        myreduce = Reduce()
        myreduce.files.extend(selected_reduced)
        myreduce.drpkg = 'maroonxdr'
        myreduce.recipename = 'applyBarycentricCorrection'
        myreduce.uparms = {'barycentricCorrection:target_name': 'HD3651'}
        myreduce.runr()

Each call writes one ``*_barycor.fits`` file into ``science_dir/``. The
output carries the ``BARYCOR`` tag.

**Verify barycor-corrected files:**

.. code-block:: python

    # List all barycor-corrected files
    print(dataselect.select_data(get_files(), tags=['BARYCOR']))


Step 9: Export Reduced Science Bundle
--------------------------------------

**Purpose**: combine the BLUE and RED barycor-corrected spectra of each
science observation into a single multi-extension FITS.

The ``exportReducedBundle`` recipe groups its inputs by the original GOA
bundle name (``ARCHNAME``), then merges each (BLUE, RED) pair into one
output. Both arms for a given observation must be passed in the **same**
``Reduce`` call - the recipe raises an error if either arm is missing for
any ``ARCHNAME``.

Select every barycor output and run the recipe once:

.. code-block:: python

    # Select all barycor-corrected files (both arms)
    selected_barycor = dataselect.select_data(
        get_files('*_barycor.fits'), tags=['PROCESSED', 'BARYCOR'])

    # Combine BLUE + RED into the final bundle
    myreduce = Reduce()
    myreduce.files.extend(selected_barycor)
    myreduce.drpkg = 'maroonxdr'
    myreduce.recipename = 'exportReducedBundle'
    myreduce.runr()

The output is a single ``<ARCHNAME>_reduced.fits`` per observation in
``science_dir/`` - for this tutorial,
``N20250717M5299_reduced.fits``. This is the science-ready product.

