.. maroonx_reduction_pipeline.rst

.. _maroonx_reduction:

****************************************
MAROON-X DRP: Using the Reduce class API
****************************************

Complete Reduction Workflow
============================

This guide describes the end-to-end data reduction workflow for MAROON-X spectroscopic observations using the DRAGONS framework. 
The pipeline processes raw data through several calibration steps to produce science-ready spectra.

Overview of Reduction Steps
---------------------------

The MAROON-X reduction pipeline has the following steps that should be executed in order:

1. **Setup and dataset** - Configure environment and extract individual arm exposures
2. **Master Flats** - Create master flat field for each arm
3. **Master Darks** - Create master dark frames for each exposure time
4. **Coefficient Darks** - Calculate dark coefficients as function of exposure time and ND filter position
5. **Synthetic Darks** - Create synthetic darks to match science exposures
6. **Wavelength Calibration** - Process Etalon and LFC calibration frames
7. **Science Reduction** - Final processing of science spectra

.. note:: MAROON-X has two separate arms (BLUE and RED) that require independent 
   calibration processing.

Prerequisites
=============

Before beginning the reduction, ensure you have:

* Raw MAROON-X FITS files organized in a science directory
* DRAGONS software with ``maroonxdr`` package installed
* All required calibration types: darks, flats, wavecal frames, and science exposures

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
    │   └── processed_flat/       # Processed flat fields
    └── test_reduction.log        # Reduction log file

Dataset
-------
In this example, we will reduce the following MX set of observations:

Here is the breakdown of the files.  All the files are included in the tutorial data
package.  They can also be downloaded from the Gemini Observatory Archive (GOA).

+---------------------+----------------------------------------------+
| Darks               || N20241115-16*                               |
+---------------------+----------------------------------------------+
| Flats               || N20241114*                                  |
+---------------------+----------------------------------------------+
| Wavecal             || N20241124*                                  |
+---------------------+----------------------------------------------+
| Science             || N20241124*                                  |
+---------------------+----------------------------------------------+




Step-by-Step Reduction Process
==============================



Step 0: Environment Setup
--------------------------

.. note:: This step will change once caldb is fully implemented.

First we import the required libraries and configure the logging file at debug level.

.. code-block:: python
    :linenos:
    :lineno-start: 1

    import os
    import itertools as it
    from pathlib import Path

    from gempy.adlibrary import dataselect
    from gempy.utils import logutils
    from recipe_system.reduction.coreReduce import Reduce

    import astrodata
    import maroonx_instruments  # noqa : important to load adclass tags

    # Configure logging
    logutils.config(file_name="test_reduction.log", mode="debug")

We then define some helper functions to get files and tags. Some steps
during the reduction will require selecting files based on exposure time
and arm. We define the tags here for easy selection.

For this particular example the exposure times are 60s, 120s, 300s, 600s, 900s, 1200s and 1800s.

.. code-block:: python
    :linenos:
    :lineno-start: 14

    # Get all files in the science_dir. Change the path here to suit your installation.
    science_dir = Path('/path/to/science_dir')
    os.chdir(science_dir)

    # Calibration files paths
    proc_dark = science_dir / "calibrations" / "processed_dark"
    proc_flat = science_dir / "calibrations" / "processed_flat"

    def get_files(path=None):
        """Get all files in a directory."""
        # If path is not provided, use the current working directory
        if path is None:
            path = science_dir
        all_files = list(path.glob('*.fits'))
        all_files = [str(f) for f in all_files]
        all_files.sort()
        return all_files

    # Define tags for file selection
    exptime_tags = ['60s', '120s', '300s', '600s', '900s', '1200s', '1800s']
    arm_tags = ['BLUE', 'RED']



Step 1: Debundle Raw Data
-------------------------

**Purpose**: Extract individual exposures from MAROON-X bundle files.

MAROON-X data is typically stored in bundle format containing RED and BLUE exposures. This step:

* Selects all files tagged as ``BUNDLE``
* Runs the automatic debundling recipe
* Creates individual FITS files for each arm color. 

This will generally produce two files per bundle, for example, ``N20241114M3295.fits`` will yield: 
``20241114T181815Z_DFFFD_b_0008.fits`` and ``20241114T181815Z_DFFFD_r_0002.fits``.

.. code-block:: python
    :linenos:
    :lineno-start: 14

    # Select bundles
    selected_bundles = dataselect.select_data(get_files(), tags=['BUNDLE'])

    # Run reduce on all selected files
    # Dragons will choose the correct recipe based on the tags
    myreduce = Reduce()
    myreduce.files.extend(selected_bundles)
    myreduce.drpkg = 'maroonxdr'
    myreduce.runr()


Step 2: Master Flat Fields
---------------------------

**Purpose**: Create master flat field calibrations.

This step processes flat field calibrations separately for each arm:

* Selects raw flat exposures (``RAW``, ``FLAT``) for BLUE and RED arms
* Processes both flat modes taken with MX: ``DFFFD`` and ``DDDDF``. This yields the combined ``DFFFF`` masterflat.

.. code-block:: python
    :linenos:
    :lineno-start: 23

    # Flats should be run for red and blue arms separately
    for arm in arm_tags:

        # Select both DFFFD and DDDDF files
        selected_flats = dataselect.select_data(get_files(), tags=['RAW', 'FLAT', arm])

        # Run reduce on all selected files
        myreduce = Reduce()
        myreduce.files.extend(selected_flats)
        myreduce.drpkg = 'maroonxdr'
        myreduce.recipename = 'makeProcessedFlatDFFFF'
        myreduce.runr()

.. note:: In the above step we use the non-default recipe **makeProcessedFlatDFFFF** to show how to process a set of flat 
    files that don't have the first fiber illuminated.

Step 3: Master Dark Frames
---------------------------

**Purpose**: Create master dark calibration frames.

This step creates master darks for all combinations of:

* **Exposure times**: 60s, 120s, 300s, 600s, 900s, 1200s, 1800s
* **Arms**: BLUE, RED

For each combination, the pipeline:

* Selects raw dark exposures matching the exposure time and arm
* Combines multiple darks to create a master dark frame
* Stores results in ``calibrations/processed_dark/``

.. code-block:: python
    :linenos:
    :lineno-start: 35

    # Dark frames should be run for each exposure time and arm
    for exptime, arm in it.product(exptime_tags, arm_tags):

        selected_darks = dataselect.select_data(get_files(), tags=['RAW', 'DARK', exptime, arm])

        # Run reduce on all selected files
        myreduce = Reduce()
        myreduce.files.extend(selected_darks)
        myreduce.drpkg = 'maroonxdr'
        myreduce.runr()

Step 4: Dark Coefficient Generation
-----------------------------------

**Purpose**: Create scaling coefficients for interpolating dark corrections between exposure times.

This step analyzes the master darks to create coefficients that allow:

* Interpolation of dark corrections as function of exposure times
* More accurate dark subtraction for science exposures

The recipe we run is called ``makeDarkCoefficients``:

* Selects all processed master darks for each arm
* Fits scaling relationships between exposure time and dark level
* Creates coefficient files for future synthetic dark generation

.. code-block:: python
    :linenos:
    :lineno-start: 35

    for arm in arm_tags:

        # Select master darks with PROCESSED tag 
        selected_darks = dataselect.select_data(
            get_files(proc_dark), tags=['PROCESSED', 'DARK', arm])

        # Run reduce on all selected files
        myreduce = Reduce()
        myreduce.files.extend(selected_darks)
        myreduce.drpkg = 'maroonxdr'
        myreduce.recipename = 'makeDarkCoefficients'
        myreduce.runr()


Step 5: Synthetic Dark Creation
-------------------------------

**Purpose**: Generate custom dark frames matched to science exposure parameters.

For science exposures, this step:

* Selects raw science exposures for each arm
* Uses dark coefficients to create synthetic darks
* Matches dark correction to exact science exposure conditions

This approach provides better dark correction than using only discrete master darks.

**Recipe**: ``makeSyntheticDark``

Step 6: Wavelength Calibration
------------------------------

**Purpose**: Establish wavelength solutions using ThAr calibration lamps.

This step processes wavelength calibration frames:

* Selects raw ThAr exposures (``WAVECAL``) for each arm
* Identifies spectral lines and fits wavelength solutions
* Creates calibration files for science spectrum wavelength assignment

**Recipe**: ``makeDynamicWavecal``

.. note:: MAROON-X uses thorium-argon (ThAr) hollow cathode lamps for wavelength 
   calibration, providing dense line coverage across the spectral range.

Step 7: Science Data Reduction
------------------------------

**Purpose**: Final processing of science spectra using all calibrations.

This step applies all calibrations to science data:

* Selects raw science exposures for each arm
* Applies dark correction, flat fielding, and wavelength calibration
* Produces calibrated, wavelength-corrected science spectra
* Example shown processes 300s exposures only

**Recipe**: Automatic recipe selection based on science tags

