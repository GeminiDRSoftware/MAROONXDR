.. maroonx_reduction_cli.rst

.. _maroonx_reduction_cli:

*********************************
MAROON-X DRP: Using the CLI tools
*********************************

Complete Reduction Workflow Using CLI
======================================

This guide walks through the end-to-end reduction of a single MAROON-X
science exposure using the DRAGONS command-line tools (``reduce``,
``dataselect``, ``caldb``). The example reduces an exposure of object
HD 3651 taken on 2025-07-17. The final output is a fully reduced,
wavelength-calibrated spectrum with science quality.

The same workflow expressed through the Python API is in
:ref:`maroonx_reduction`. Both tutorials reduce the same files and produce
the same outputs; they differ only in interface.

Overview of Reduction Steps
----------------------------

The full reduction proceeds in nine steps, executed in order:

1. **Debundle** - split raw GOA bundles into per-arm FITS files
2. **Master Darks** - combine raw darks per exposure time and arm
3. **Dark Coefficients** - fit dark scaling vs. exposure time per arm
4. **Master Flats** - build per-arm master flats from DFFFD + DDDDF inputs
5. **Wavelength Calibration** - fit the dynamic etalon wavelength solution
6. **Synthetic Darks** - interpolate a dark matched to each science exposure
7. **Science Reduction** - extract and wavelength-calibrate science spectra
8. **Barycentric Correction** - apply the BERV shift to the reduced spectra
9. **Export Bundle** - combine BLUE and RED arms into a single output FITS

.. note:: MAROON-X has two independent arms (BLUE and RED). Every
   calibration and science step is run separately per arm.

.. important:: MAROON-X is currently an external instrument package, so all
   DRAGONS CLI tools require ``--adpkg maroonx_instruments`` (and
   ``--drpkg maroonxdr`` for ``reduce``) to recognise the file tags and
   load the recipes. These flags will no longer be needed once MAROON-X is
   integrated into DRAGONS.

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

As the reduction proceeds, ``reduce`` writes its outputs into the working
directory and registers calibrations under a ``calibrations/`` subtree
managed by ``caldb``::

    science_dir/
    ├── N*.fits                          # Raw GOA bundle files
    ├── *_b_*.fits, *_r_*.fits           # Debundled per-arm files and
    │                                    # downstream reduction products
    ├── calibrations/
    │   ├── processed_dark/              # Master darks (per exptime, per arm)
    │   ├── processed_dark_coeff/        # Dark scaling coefficients (per arm)
    │   ├── processed_flat/              # Master flats (per arm)
    │   └── processed_arc/               # Dynamic etalon wavelength solutions
    └── reduce_*.log                     # Reduction log files

.. note:: Calibration files produced by ``reduce`` are written twice: once
   in ``science_dir/`` and once in the matching ``calibrations/<caltype>/``
   subdirectory. The copy under ``calibrations/`` is what ``caldb`` indexes
   and serves to later steps.

.. todo:: Link to the DRAGONS ``caldb`` setup page (initialising the
   calibration database, ``~/.dragonsrc`` configuration, etc.) once that
   section of the User Manual is written.

Dataset
-------

The worked example reduces a single 300 s exposure of HD 3651 in the SOOOE
fiber pattern, taken on 2025-07-17, against calibrations taken on
2025-07-01, 2025-07-17, and 2025-07-21.

All files can be downloaded from the Gemini Observatory Archive (GOA). The
table below lists the raw GOA bundle ranges; each ``N…M*.fits`` bundle
contains one BLUE and one RED exposure that ``splitBundle`` (Step 1)
separates into individual files.

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
=================


Step 1: Debundle Raw Data
--------------------------

**Purpose**: split each raw GOA bundle into one BLUE and one RED FITS file.

MAROON-X data arrives from GOA as bundle files (tagged ``BUNDLE``) that
contain both arms in a single FITS. The default ``processBundle`` recipe
picks up every ``BUNDLE``-tagged file and writes a per-arm FITS for each
one. For example, the science bundle ``N20250717M5299.fits`` yields
``20250717T144308Z_SOOOE_b_0300.fits`` and
``20250717T144308Z_SOOOE_r_0004.fits``.

Select every bundle in the working directory and run ``reduce``:

.. code-block:: bash

    # Build the bundle list
    dataselect --adpkg maroonx_instruments --tags BUNDLE -o bundles.lis *.fits

    # Debundle (no --recipe needed: processBundle is the default for BUNDLE)
    reduce --adpkg maroonx_instruments --drpkg maroonxdr @bundles.lis

.. note:: The ``@`` prefix tells ``reduce`` to read filenames from the
   given list file, one per line.

After this step ``science_dir/`` contains the original ``N*.fits`` bundles
plus the per-arm debundled files used by every subsequent step.


Step 2: Master Darks
---------------------

**Purpose**: combine the raw darks at each exposure time into a master
dark, per arm.

For every (exposure time, arm) pair, ``dataselect`` finds the matching
raw darks (DDDDE pattern, tagged ``RAW,DARK``) and ``reduce`` runs the
default ``makeProcessedDark`` recipe. The dataset has 7 exposure times ×
2 arms, so 14 master darks total.

For a single combination of exposure time and arm, the process is:

.. code-block:: bash

    # Process darks for a specific exposure time and arm, 300s and BLUE
    dataselect --adpkg maroonx_instruments --tags RAW,DARK,300s,BLUE \
        -o darks_300s_BLUE.lis *.fits
    reduce --adpkg maroonx_instruments --drpkg maroonxdr @darks_300s_BLUE.lis

Run the full loop directly in the shell:

.. code-block:: bash

    # Process darks for each exposure time and arm
    for exptime in 60s 120s 300s 600s 900s 1200s 1800s; do
        for arm in BLUE RED; do

            # Select darks for this combination and save to list
            dataselect --adpkg maroonx_instruments --tags RAW,DARK,$exptime,$arm \
                -o darks_${exptime}_${arm}.lis *.fits

            # Reduce the list
            reduce --adpkg maroonx_instruments --drpkg maroonxdr \
                @darks_${exptime}_${arm}.lis
        done
    done

Each output is a single master dark file named after the first raw frame
in the stack with a ``_dark`` suffix - for example
``20250721T170049Z_DDDDE_b_0300_dark.fits``. The 14 master darks land in
``calibrations/processed_dark/`` (and as duplicates in ``science_dir/``).

**Verify processed darks:**

.. code-block:: bash

    # List all processed darks
    dataselect --adpkg maroonx_instruments --tags PROCESSED,DARK *.fits

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
input selector takes all ``PROCESSED,DARK`` files for the arm and excludes
the ``DARK_COEFF`` and ``DARK_SYNTH`` tags so previous coefficient or
synthetic-dark outputs in the same directory are not re-ingested.

Run once per arm:

.. code-block:: bash

    for arm in BLUE RED; do

        # Select all master darks for this arm, excluding any previous
        # coefficient or synthetic-dark outputs
        dataselect --adpkg maroonx_instruments \
            --tags PROCESSED,DARK,$arm \
            --xtags DARK_COEFF,DARK_SYNTH \
            -o dark_coeffs_${arm}.lis *.fits

        # Fit the dark-coefficient model
        reduce --adpkg maroonx_instruments --drpkg maroonxdr \
            --recipe makeDarkCoefficients @dark_coeffs_${arm}.lis
    done

Each call writes one ``*_darkCoefficients.fits`` file (named after the
first input master dark) into ``science_dir/``, and a copy under the
``processed_dark_coeff`` caltype in ``calibrations/processed_dark_coeff/``.
The output carries the ``DARK_COEFF`` tag, which is how later steps locate
it.

**Verify dark-coefficient files:**

.. code-block:: bash

    # Select dark coefficient files
    dataselect --adpkg maroonx_instruments --tags DARK_COEFF *.fits


Step 4: Master Flats
---------------------

**Purpose**: build a master flat per arm.

MAROON-X flats come in two illumination patterns: ``DFFFD`` (flat on
fibers 2-4) and ``DDDDF`` (flat on fiber 5 only). The dataset for this
tutorial contains both patterns and uses the ``makeProcessedFlatDFFFF``
recipe, which combines the two into a single ``DFFFF`` master flat per
arm.

Run the reduction per arm:

.. code-block:: bash

    for arm in BLUE RED; do

        # Select every raw flat for this arm (both DFFFD and DDDDF)
        dataselect --adpkg maroonx_instruments --tags RAW,FLAT,$arm \
            -o flats_${arm}.lis *.fits

        # Build the master flat
        reduce --adpkg maroonx_instruments --drpkg maroonxdr \
            --recipe makeProcessedFlatDFFFF @flats_${arm}.lis
    done

Each call writes one ``*_flat.fits`` master flat into ``science_dir/`` and
a copy under the ``processed_flat`` caltype in
``calibrations/processed_flat/``. The output carries the ``PROCESSED,FLAT``
tag set.

.. note:: The default recipe ``makeProcessedFlat`` combines ``DFFFD`` and
   ``FDDDF`` patterns into an ``FFFFF`` master flat. When the dataset has
   ``DDDDF`` instead of ``FDDDF`` (as in this tutorial), use the
   ``makeProcessedFlatDFFFF`` variant explicitly via ``--recipe``.

**Verify processed flats:**

.. code-block:: bash

    # List all processed flats
    dataselect --adpkg maroonx_instruments --tags PROCESSED,FLAT *.fits


Step 5: Wavelength Calibration
-------------------------------

**Purpose**: fit, per arm, a dynamic wavelength solution from the etalon
calibration frame.

The Fabry-Perot etalon is the primary wavelength calibrator for MAROON-X.
Anchored to an initial ThAr solution, the etalon provides the refined
wavelength and drift solution that the science reduction (Step 7) applies
to extracted spectra.

The ``makeDynamicWavecal`` recipe takes the raw etalon frame
(``DEEEE`` pattern, tagged ``RAW,ETALON``) for one arm, identifies the
peaks in every order, and fits the per-pixel wavelength solution. The
``getPeaksAndPolynomials:multithreading=True`` override parallelizes the
fit and is recommended on multi-core machines.

Run once per arm:

.. code-block:: bash

    for arm in BLUE RED; do

        # Select the raw etalon wavecal frame for this arm
        dataselect --adpkg maroonx_instruments --tags RAW,ETALON,$arm \
            -o wavecal_${arm}.lis *.fits

        # Fit the dynamic wavelength solution
        reduce --adpkg maroonx_instruments --drpkg maroonxdr \
            --recipe makeDynamicWavecal \
            -p getPeaksAndPolynomials:multithreading=True \
            @wavecal_${arm}.lis
    done

Each call writes one ``*_wavecal.fits`` file into ``science_dir/`` and a
copy under the ``processed_arc`` caltype in ``calibrations/processed_arc/``.
The output carries the ``PROCESSED,ARC`` tag set.

**Verify wavelength calibrations:**

.. code-block:: bash

    # List processed wavecal files
    dataselect --adpkg maroonx_instruments --tags PROCESSED,ARC *.fits


Step 6: Synthetic Darks
------------------------

**Purpose**: interpolate, per science frame, a synthetic dark matched to
its exact exposure time using the dark coefficients fit in Step 3.

The ``makeSyntheticDark`` recipe takes the raw science frames and uses the
``DARK_COEFF`` calibration to build a dark frame at the science exposure
time. The synthetic darks are stored under the ``processed_dark`` caltype
with the ``DARK_SYNTH`` tag so Step 7 can find them as ordinary darks.

Run once per arm:

.. code-block:: bash

    for arm in BLUE RED; do

        # Select raw science frames for this arm
        dataselect --adpkg maroonx_instruments --tags RAW,SCI,$arm \
            -o sci_${arm}.lis *.fits

        # Build the synthetic dark for each science frame
        reduce --adpkg maroonx_instruments --drpkg maroonxdr \
            --recipe makeSyntheticDark @sci_${arm}.lis
    done


Step 7: Science Data Reduction
-------------------------------

**Purpose**: extract and wavelength-calibrate the science spectra, applying
the master flat, the dynamic wavelength solution, and the synthetic dark
matched to the science exposure time.

The default recipe for ``RAW,SCI`` files runs the full echelle reduction:
dark and flat correction, straylight removal, stripe extraction, optimal
extraction, wavelength assignment, and fiber combination. ``caldb``
resolves each input calibration automatically from the work done in
Steps 2-6.

The three ``-p`` overrides below are the recommended defaults for the
current pipeline: ``extractStripes:straylight_removal_fibers=[5]`` applies
the straylight correction to the calibration fiber, ``getPeaksAndPolynomials``
runs the per-order fits in parallel, and ``combineFibers:max_clips=20``
sets the outlier rejection threshold when combining the science fibers.

Run once per arm:

.. code-block:: bash

    for arm in BLUE RED; do

        # Select raw science frames for this arm
        dataselect --adpkg maroonx_instruments --tags RAW,SCI,$arm \
            -o sci_${arm}.lis *.fits

        # Full echelle reduction
        reduce --adpkg maroonx_instruments --drpkg maroonxdr \
            -p 'extractStripes:straylight_removal_fibers=[5]' \
            -p 'getPeaksAndPolynomials:multithreading=True' \
            -p 'combineFibers:max_clips=20' \
            @sci_${arm}.lis
    done

Each call writes a ``*_reduced.fits`` file per input science frame into
``science_dir/``. The reduced products carry the ``PROCESSED_SCIENCE`` tag
and are the input to the barycentric-correction step (Step 8).

**Verify reduced science frames:**

.. code-block:: bash

    # List all reduced science frames
    dataselect --adpkg maroonx_instruments --tags PROCESSED_SCIENCE *.fits


Step 8: Barycentric Correction
-------------------------------

**Purpose**: shift the reduced spectra into the solar-system barycentric
rest frame.

The ``applyBarycentricCorrection`` recipe reads the time, sky-position,
and observatory descriptors from each reduced file, queries SIMBAD to
resolve the target's coordinates, computes the barycentric Earth radial
velocity (BERV), and writes the corrected spectrum to ``*_barycor.fits``.

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

.. code-block:: bash

    for arm in BLUE RED; do

        # Select the reduced science file for this arm
        dataselect --adpkg maroonx_instruments --tags PROCESSED_SCIENCE,$arm \
            -o barycor_${arm}.lis *_reduced.fits

        # Apply the barycentric correction
        reduce --adpkg maroonx_instruments --drpkg maroonxdr \
            --recipe applyBarycentricCorrection \
            -p barycentricCorrection:target_name=HD3651 \
            @barycor_${arm}.lis
    done

Each call writes one ``*_barycor.fits`` file into ``science_dir/``. The
output carries the ``BARYCOR`` and ``PROCESSED_SCIENCE`` tags.

**Verify barycor-corrected files:**

.. code-block:: bash

    # List all barycor-corrected files
    dataselect --adpkg maroonx_instruments --tags PROCESSED,BARYCOR *.fits


Step 9: Export Reduced Science Bundle
--------------------------------------

**Purpose**: combine the BLUE and RED barycor-corrected spectra of each
science observation into a single multi-extension FITS.

The ``exportReducedBundle`` recipe groups its inputs by the original GOA
bundle name (``ARCHNAME``), then merges each (BLUE, RED) pair into one
output. Both arms for a given observation must be passed in the **same**
``reduce`` call - the recipe raises an error if either arm is missing for
any ``ARCHNAME``.

Select every barycor output and run the recipe once:

.. code-block:: bash

    # Select all barycor-corrected files (both arms)
    dataselect --adpkg maroonx_instruments --tags PROCESSED,BARYCOR \
        -o bundle.lis *_barycor.fits

    # Combine BLUE + RED into the final bundle
    reduce --adpkg maroonx_instruments --drpkg maroonxdr \
        --recipe exportReducedBundle @bundle.lis

The output is a single ``<ARCHNAME>_reduced.fits`` per observation in
``science_dir/`` - for this tutorial,
``N20250717M5299_reduced.fits``. This is the science-ready product.


Advanced CLI Usage
==================

This section covers a few CLI tools and patterns that go beyond the
linear walk-through above. They are useful when something goes wrong, or
when you want to inspect the pipeline state interactively.

Inspecting files: ``typewalk``
-------------------------------

``typewalk`` prints the AstroData tag set of one or more FITS files. It is
the first thing to reach for when ``dataselect`` is not picking up the
files you expect:

.. code-block:: bash

    # Show all tags on a single file
    typewalk --adpkg maroonx_instruments N20250717M5299.fits

    # Walk a whole directory and print tags per file
    typewalk --adpkg maroonx_instruments .

Tag combinations such as ``RAW``, ``SCI``, ``BLUE``, ``300s``,
``PROCESSED_SCIENCE``, ``BARYCOR`` are exactly the strings ``dataselect``
will match against.

Discovering recipes and parameters: ``showrecipes`` / ``showpars``
-------------------------------------------------------------------

The DRAGONS framework provides two discovery commands that answer "what
recipes can run on this file?" and "what parameters does this primitive
take?":

.. code-block:: bash

    # List every recipe applicable to this file's tag set
    showrecipes --adpkg maroonx_instruments --drpkg maroonxdr \
        20250717T144308Z_SOOOE_b_0300.fits

    # Show the default parameters for a primitive (and its docstring)
    showpars --adpkg maroonx_instruments --drpkg maroonxdr \
        20250717T144308Z_SOOOE_b_0300.fits combineFibers

Use these to discover the recipe names you can pass to ``reduce --recipe``
and the parameter names you can pass to ``-p primitive:param=value``.

Overriding primitive parameters
--------------------------------

The ``-p primitive:param=value`` syntax used throughout this tutorial sets
a single parameter on a single primitive for the current ``reduce`` call.
Several ``-p`` flags can be chained:

.. code-block:: bash

    reduce --adpkg maroonx_instruments --drpkg maroonxdr \
        -p extractStripes:straylight_removal_fibers=[5] \
        -p getPeaksAndPolynomials:multithreading=True \
        -p combineFibers:max_clips=20 \
        @sci_BLUE.lis

Discover available parameters with ``showpars`` (above).

Managing the calibration database: ``caldb``
---------------------------------------------

``reduce`` registers processed calibrations with ``caldb`` automatically,
so the linear walk-through above never needs to call ``caldb`` directly.
For troubleshooting - "this run picked the wrong flat", "the wavecal from
yesterday is shadowing today's" - ``caldb`` has a small set of commands:

.. code-block:: bash

    caldb list                            # show every calibration registered
    caldb add calibrations/processed_flat/<file>.fits
    caldb remove <file>.fits

.. todo:: Link to the full ``caldb`` documentation page once that section
   of the User Manual is written.

Logging
-------

Two flags control ``reduce`` logging:

.. code-block:: bash

    # Set log verbosity (mode: quiet | standard | debug)
    reduce --adpkg maroonx_instruments --drpkg maroonxdr \
        --logmode debug @files.lis

    # Send the log to a specific file (default: reduce.log)
    reduce --adpkg maroonx_instruments --drpkg maroonxdr \
        --logfile my_reduction.log @files.lis

