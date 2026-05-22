.. maroonx_reduction_cli.rst

.. _maroonx_reduction_cli:

*********************************
MAROON-X DRP: Using the CLI tools
*********************************

Complete Reduction Workflow Using CLI
======================================

This guide describes the end-to-end data reduction workflow for MAROON-X spectroscopic observations using DRAGONS command-line tools.
The pipeline processes raw data through several calibration steps to produce science-ready spectra.

This tutorial demonstrates the same reduction workflow as :ref:`maroonx_reduction` but uses command-line tools instead of the Python API.

Overview of Reduction Steps
----------------------------

The MAROON-X reduction pipeline has the following steps that should be executed in order:

1. **Setup and dataset** - Configure environment and extract individual arm exposures
2. **Master Flats** - Create master flat field for each arm
3. **Master Darks** - Create master dark frames for each exposure time
4. **Coefficient Darks** - Calculate dark coefficients as function of exposure time and ND filter position
5. **Synthetic Darks** - Create synthetic darks to match science exposures
6. **Wavelength Calibration** - Process Etalon/LFC calibration frames
7. **Science Reduction** - Final processing of science spectra

.. note:: MAROON-X has two separate arms (BLUE and RED) that require independent
   calibration processing.

.. important:: Since MAROON-X is developed as an external instrument package, all DRAGONS CLI tools
   require the ``--adpkg maroonx_instruments`` flag to properly recognize MAROON-X files and tags.
   This should not be required when MAROON-X is fully integrated into DRAGONS.

Prerequisites
=============

Before beginning the reduction, ensure you have:

* Raw MAROON-X FITS files organized in a science directory
* DRAGONS software with ``maroonxdr`` package installed, preferably in a virtual environment.

Directory Structure
-------------------

The example assumes the following directory structure::

    science_dir/
    └── *.fits                    # Raw science and calibration files

As the reduction proceeds, the following subdirectories will be created::

    science_dir/
    ├── *.fits                    # Raw science and calibration files
    ├── calibrations/
    │   ├── processed_dark/       # Processed dark frames
    │   ├── processed_flat/       # Processed flat fields
    │   └── processed_wavecal/    # Processed wavecal frames
    └── reduce_*.log              # Reduction log files

Dataset
-------
In this example, we will reduce the following MX set of observations:

Here is the breakdown of the files. They can be downloaded from the Gemini Observatory Archive (GOA).

+---------------------+----------------------------------------------+
| Darks               || N20241115-16*                               |
+---------------------+----------------------------------------------+
| Flats               || N20241114*                                  |
+---------------------+----------------------------------------------+
| Wavecal             || N20241124*                                  |
+---------------------+----------------------------------------------+
| Science             || N20241124*                                  |
+---------------------+----------------------------------------------+


For this particular example, the exposure times are: **60s, 120s, 300s, 600s, 900s, 1200s, and 1800s**.

.. todo:: Describe the dataset and their tags used to select them.


Reduction Steps
=================


Step 1: Debundle Raw Data
--------------------------

**Purpose**: Extract individual exposures from MAROON-X bundle files.

MAROON-X data is typically stored in bundle format containing RED and BLUE exposures. This step:

* Selects all files tagged as ``BUNDLE``
* Runs the automatic debundling recipe
* Creates individual FITS files for each arm color

This will generally produce two files per bundle, for example, ``N20241114M3295.fits`` will yield:
``20241114T181815Z_DFFFD_b_0008.fits`` and ``20241114T181815Z_DFFFD_r_0002.fits``.

**Using dataselect to identify bundle files:**

.. code-block:: bash

    # List all bundle files
    dataselect --adpkg maroonx_instruments --tags BUNDLE *.fits

**Run the reduction:**

.. code-block:: bash

    # Create file list
    dataselect --adpkg maroonx_instruments --tags BUNDLE -o bundles.lis *.fits

    # Run reduce on all bundles
    reduce --adpkg maroonx_instruments --drpkg maroonxdr @bundles.lis

.. note:: The ``@`` symbol tells reduce to read filenames from the list file.


Step 2: Master Flat Fields
---------------------------

**Purpose**: Create master flat field calibrations.

This step processes flat field calibrations separately for each arm:

* Selects raw flat exposures: ``RAW``, ``FLAT``, ``BLUE`` (or ``RED``).
* Processes both flat modes taken with MX: ``DFFFD`` and ``FDDDF`` . This yields the combined ``FFFFF`` masterflat.
However due to problems illuminating the first fiber, we could have ``DDDDF`` instead of ``FDDDF``.
For this case we use the non-default recipe **makeProcessedFlatDFFFF** (``--recipe``, ``-r``).

**Process BLUE arm flats:**

.. code-block:: bash

    # Select BLUE arm flats
    dataselect --adpkg maroonx_instruments --tags RAW,FLAT,BLUE -o flats_blue.lis *.fits

    # Reduce with default recipe: DFFFD + FDDDF -> FFFFF
    reduce --adpkg maroonx_instruments --drpkg maroonxdr @flats_blue.lis
    
    # Reduce with specific recipe: DFFFD + DDDDF -> DFFFF
    reduce --adpkg maroonx_instruments --drpkg maroonxdr --recipe makeProcessedFlatDFFFF @flats_blue.lis


**Process RED arm flats:**

.. code-block:: bash

    # Select RED arm flats
    dataselect --adpkg maroonx_instruments --tags RAW,FLAT,RED -o flats_red.lis *.fits

    # Reduce with default recipe: DFFFD + FDDDF -> FFFFF
    reduce --adpkg maroonx_instruments --drpkg maroonxdr @flats_red.lis

    # Reduce with specific recipe: DFFFD + DDDDF -> DFFFF
    reduce --adpkg maroonx_instruments --drpkg maroonxdr --recipe makeProcessedFlatDFFFF @flats_red.lis

**Verify processed flats:**

After reduction is complete the master flat files should be present in the calibration directory.

.. code-block:: bash

    # List all processed flats
    dataselect --adpkg maroonx_instruments --tags PROCESSED,FLAT calibrations/processed_flat/*.fits

.. note::
    Processed calibration files will be saved by duplicate in the current working directory ``science_dir/``
    and in the respective calibration directories.

Step 3: Master Dark Frames
---------------------------

**Purpose**: Create master dark calibration frames.

This step creates master darks using 2 to 5 DDDDE frames per exposure time and per arm, in this case:

* **Exposure times**: 60s, 120s, 300s, 600s, 900s, 1200s, 1800s
* **Arms**: BLUE, RED

For a single combination of exposure time and arm, the process is:

.. code-block:: bash

    # Process darks for a specific exposure time and arm, 300s and BLUE
    dataselect --adpkg maroonx_instruments --tags RAW,DARK,300s,BLUE -o darks_300s_blue.lis *.fits
    reduce --adpkg maroonx_instruments --drpkg maroonxdr @darks_300s_blue.lis

**Process all dark combinations:**

.. code-block:: bash

    # Process darks for each exposure time and arm
    for exptime in 60s 120s 300s 600s 900s 1200s 1800s; do
        for arm in BLUE RED; do

            # Select darks for this combination
            darks=$(dataselect --adpkg maroonx_instruments --tags RAW,DARK,$exptime,$arm *.fits)

            # Reduce the selected darks
            reduce --adpkg maroonx_instruments --drpkg maroonxdr $darks
        done
    done

**Alternative - using file lists:**

.. code-block:: bash

    # Create separate list files for each combination
    for exptime in 60s 120s 300s 600s 900s 1200s 1800s; do
        for arm in BLUE RED; do

            # Select darks for this combination and save to list
            dataselect --adpkg maroonx_instruments --tags RAW,DARK,$exptime,$arm -o darks_${exptime}_${arm}.lis *.fits

            # Reduce the list
            reduce --adpkg maroonx_instruments --drpkg maroonxdr @darks_${exptime}_${arm}.lis
        done
    done

**Verify processed darks:**

After reduction is complete the master dark files should be present in the calibration directory.

.. code-block:: bash

    # List all processed darks
    dataselect --adpkg maroonx_instruments --tags PROCESSED,DARK calibrations/processed_dark/*.fits

.. note::
    Processed calibration files will be saved by duplicate in the current working directory ``science_dir/``
    and in the respective calibration directories.

Step 4: Dark Coefficient Generation
------------------------------------

**Purpose**: Create scaling coefficients for interpolating dark corrections between exposure times.

This step analyzes the master darks to create coefficients that allow:

* Interpolation of dark corrections as function of exposure times
* Creation of synthetic dark for specific science exposures

**Process BLUE arm coefficients:**

The recipe we run is called ``makeDarkCoefficients``:

.. code-block:: bash

    # Select all processed BLUE darks
    dataselect --adpkg maroonx_instruments --tags PROCESSED,DARK,BLUE \
        -o proc_darks_blue.lis calibrations/processed_dark/*.fits

    # Generate coefficients
    reduce --adpkg maroonx_instruments --drpkg maroonxdr \
        -r makeDarkCoefficients @proc_darks_blue.lis

**Process RED arm coefficients:**

.. code-block:: bash

    # Select all processed RED darks
    dataselect --adpkg maroonx_instruments --tags PROCESSED,DARK,RED \
        -o proc_darks_red.lis calibrations/processed_dark/*.fits

    # Generate coefficients
    reduce --adpkg maroonx_instruments --drpkg maroonxdr \
        -r makeDarkCoefficients @proc_darks_red.lis

This reduction will generate new FITS files that are stored and interpreted as processed darks, and are 
saved in ``calibrations/processed_dark``. To make them detectable they have the tag ``DARK_COEFF``.

Verify new coefficient files using the tag ``DARK_COEFF``:

.. code-block:: bash

    # Select dark coefficient files
    dataselect --adpkg maroonx_instruments --tags DARK_COEFF calibrations/processed_dark/*.fits



Step 5: Synthetic Dark Creation
--------------------------------

**Purpose**: Generate custom dark frames matched to science exposure parameters.

For this step we select the raw science frames that we want to reduce and pass them
to the recipe ``makeSyntheticDark``. The interpolation is done using the Dark Coefficient
files created previously.

**Create BLUE arm synthetic darks:**

.. code-block:: bash

    # Select blue science frames
    dataselect --adpkg maroonx_instruments --tags RAW,SCI,BLUE -o sci_blue.lis *.fits

    # Reduce with makeSyntheticDark
    reduce --adpkg maroonx_instruments --drpkg maroonxdr -r makeSyntheticDark @sci_blue.lis

**Create RED arm synthetic darks:**

.. code-block:: bash

    # Select red science frames
    dataselect --adpkg maroonx_instruments --tags RAW,SCI,RED -o sci_red.lis *.fits

    # Reduce with makeSyntheticDark
    reduce --adpkg maroonx_instruments --drpkg maroonxdr -r makeSyntheticDark @sci_red.lis


Step 6: Wavelength Calibration
-------------------------------

**Purpose**: Establish wavelength solutions using etalon calibration frames.

The Fabry-Perot Etalon is the primary wavelength calibrator for MAROON-X. 
Based on an initial solution anchored to ThAr, the etalon provides a refined wavelength and 
drift solution for the instrument. Wavelength and drift solution follow a multistep process.

* Selects raw wavelength calibration exposures (``WAVECAL``) for each arm
* Identifies spectral peaks and fits wavelength solutions
* Creates calibration files for science spectrum wavelength assignment


**Process BLUE arm wavelength calibration:**

.. code-block:: bash

    # Select BLUE wavecal frames
    dataselect --adpkg maroonx_instruments --tags RAW,WAVECAL,BLUE -o wavecal_blue.lis *.fits

    # Generate wavelength solution
    reduce --adpkg maroonx_instruments --drpkg maroonxdr @wavecal_blue.lis

**Process RED arm wavelength calibration:**

.. code-block:: bash

    # Select RED wavecal frames
    dataselect --adpkg maroonx_instruments --tags RAW,WAVECAL,RED -o wavecal_red.lis *.fits

    # Generate wavelength solution
    reduce --adpkg maroonx_instruments --drpkg maroonxdr @wavecal_red.lis

**Verify wavelength calibrations:**

These calibrations should be stored in the calibration database (and in the cwd as duplicates).

.. code-block:: bash

    # List processed wavecal files
    dataselect --adpkg maroonx_instruments --tags PROCESSED,WAVECAL *.fits

Step 7: Science Data Reduction
-------------------------------

**Purpose**: Final processing of science spectra using all calibrations.

This step applies all calibrations to science data:

* Selects raw science exposures for each arm
* Applies dark, flat and wavelength calibrations
* Produces calibrated, wavelength-corrected science spectra

**Process BLUE arm science frames:**

.. code-block:: bash

    # Select BLUE science frames
    dataselect --adpkg maroonx_instruments --tags RAW,SCI,BLUE -o sci_blue.lis *.fits

    # Reduce science frames
    reduce --adpkg maroonx_instruments --drpkg maroonxdr @sci_blue.lis

**Process RED arm science frames:**

.. code-block:: bash

    # Select RED science frames
    dataselect --adpkg maroonx_instruments --tags RAW,SCI,RED -o sci_red.lis *.fits

    # Reduce science frames
    reduce --adpkg maroonx_instruments --drpkg maroonxdr @sci_red.lis

**Process specific exposure times:**

If you want to reduce specific exposure time set of files you can use the time exposure tag:

.. code-block:: bash

    # Select 300s, RED science frames
    dataselect --adpkg maroonx_instruments --tags RAW,SCI,RED,300s -o sci_300s_red.lis *.fits

    # Process only 300s exposures for RED arm
    reduce --adpkg maroonx_instruments --drpkg maroonxdr @sci_300s_red.lis

**Verify reduced science frames:**

.. code-block:: bash

    # List all processed science frames
    dataselect --adpkg maroonx_instruments --tags PROCESSED,SCI *.fits


Step 8: Export Reduced Science Bundles
--------------------------------------

**Purpose**: Export in a single FITS both reduced arms spectra.

.. code-block:: bash

    # Select all science frames, except bundles
    # use --xtags BUNDLE to drop previous run outputs
    dataselect --adpkg maroonx_instruments --tags PROCESSED,SCI -o sci_arms.lis *.fits

    # Reduce science frames
    reduce --adpkg maroonx_instruments --drpkg maroonxdr \
        -r exportReducedBundle @sci_arms.lis


Advanced CLI Usage
==================

Useful dataselect Patterns
---------------------------

**Select by multiple tags (AND logic):**

.. code-block:: bash

    # Files that are both DARK and BLUE and RAW
    dataselect --adpkg maroonx_instruments --tags RAW,DARK,BLUE *.fits

**Select with tag exclusions:**

.. code-block:: bash

    # All science data except bundles
    dataselect --adpkg maroonx_instruments --tags SCI --xtags BUNDLE *.fits


Reduce Command Options
-----------------------

**Specify recipe explicitly:**

.. code-block:: bash

    reduce --adpkg maroonx_instruments --drpkg maroonxdr \
        --recipe exportReducedBundle @sci_arms.lis

**Override primitive parameters:**

.. code-block:: bash

    # Override a single parameter
    reduce --adpkg maroonx_instruments --drpkg maroonxdr \
        -r exportReducedBundle \
        -p bundleArmStreams:suffix=_exported @sci_arms.lis


**Logging control:**

.. code-block:: bash

    # Set log level to debug
    reduce --adpkg maroonx_instruments --drpkg maroonxdr --logmode debug @files.lis

    # Specify log file name
    reduce --adpkg maroonx_instruments --drpkg maroonxdr --logfile my_reduction.log @files.lis

