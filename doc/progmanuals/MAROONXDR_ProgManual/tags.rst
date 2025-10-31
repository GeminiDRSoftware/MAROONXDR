.. _tags:

****
Tags
****

Tags are metadata labels automatically assigned to MAROONX data files. They determine which recipes process the data and enable data selection for calibration matching.

MAROONX uses a **fiber-configuration-based tagging system** where tags are assigned based on the 5-fiber setup pattern rather than OBSTYPE keywords.


Fiber Configuration Patterns
=============================

Tags are determined by fiber setup (5 fibers). Common patterns:

.. list-table::
   :header-rows: 1
   :widths: 12 12 12 12 12 40

   * - Fiber 1
     - Fiber 2
     - Fiber 3
     - Fiber 4
     - Fiber 5
     - Tags
   * - Flat lamp
     - Dark
     - Dark
     - Dark
     - Flat lamp
     - ``FLAT, CAL``
   * - Dark
     - Flat lamp
     - Flat lamp
     - Flat lamp
     - Dark
     - ``FLAT, CAL``
   * - Dark
     - Dark
     - Dark
     - Dark
     - Etalon
     - ``DARK, CAL``
   * - Sky
     - Target
     - Target
     - Target
     - Etalon
     - ``SCI, SPECT``
   * - Dark
     - ThAr
     - ThAr
     - ThAr
     - ThAr
     - ``WAVECAL, SPECT, ThAr, CAL``
   * - Dark
     - Etalon
     - Etalon
     - Etalon
     - Etalon
     - ``WAVECAL, SPECT, ETALON, CAL``
   * - Dark
     - LFC
     - LFC
     - LFC
     - LFC
     - ``WAVECAL, SPECT, LFC, CAL``

**Fiber types:** ``Dark``, ``Flat lamp``, ``Sky``, ``Target``, ``Etalon``, ``ThAr``, ``LFC``, ``Iodine cell``




Tag to Recipe Mapping
======================

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - Tags
     - Recipe
   * - ``MAROONX, BUNDLE``
     - :meth:`~maroonxdr.maroonx.recipes.sq.recipes_BUNDLE.processBundle`
   * - ``MAROONX, CAL, DARK``
     - :meth:`~maroonxdr.maroonx.recipes.sq.recipes_DARK.makeProcessedDark`
   * - ``MAROONX, CAL, FLAT``
     - :meth:`~maroonxdr.maroonx.recipes.sq.recipes_FLAT_SPECT.makeProcessedFlat` , :meth:`~maroonxdr.maroonx.recipes.sq.recipes_FLAT_SPECT.makeProcessedFlatDFFFF`
   * - ``MAROONX, WAVECAL, ThAr``
     - :meth:`~maroonxdr.maroonx.recipes.sq.recipes_STATICWAVECAL.makeStaticWavecal`
   * - ``MAROONX, WAVECAL``
     - :meth:`~maroonxdr.maroonx.recipes.sq.recipes_DYNAMIC_WAVECAL.makeDynamicWavecal`
   * - ``MAROONX, SCI``
     - :meth:`~maroonxdr.maroonx.recipes.sq.recipes_ECHELLE_SPECT.reduce`



Using Tags for Data Selection
==============================

List tags on files::

    typewalk *.fits

Select by tags::

    # All wavelength calibrations
    dataselect --tags WAVECAL *.fits

    # Blue arm science
    dataselect --tags SCI,BLUE *.fits

    # Processed flats
    dataselect --tags FLAT,PROCESSED *.fits

    # Raw darks with 300s exposure
    dataselect --tags RAW,DARK,300s *.fits

Exclude tags with ``--xtags``::

    # Regular darks, excluding coefficients and synthetic
    dataselect --tags DARK --xtags DARK_COEFF,DARK_SYNTH *.fits

In Python scripts::

    from recipe_system.utils import dataselect

    # Select blue arm science with 300s exposure
    sci_files = dataselect.select_data(
        all_files,
        tags=['BLUE', 'SCI', '300s']
    )


All Tags
=========

**Instrument/Observatory**

``GEMINI``, ``MAROONX``


**Arm Configuration** (mutually exclusive)

``BUNDLE`` : Raw GOA archive, blue and red arms as extensions

``BLUE`` : Blue arm data (450-670 nm)

``RED`` : Red arm data (650-920 nm)


**Observation Type**

``SPECT`` : Spectroscopic data

``SCI`` : Science observation

``CAL`` : Calibration frame

``WAVECAL`` : Calibration frame


**Calibration Types**

``DARK`` : Dark frame calibration

``FLAT`` : Flat field calibration

``WAVECAL`` : Wavelength calibration

``ThAr`` : Thorium-Argon arc lamp calibration

``ETALON`` : Etalon calibration

``LFC`` : Laser frequency comb calibration

``BPM`` : Bad pixel mask


**Processing Status**

``RAW`` : Unprocessed data

``PROCESSED`` : Processed by DRAGONS

``PREPARED`` : Basic preparation applied

**Special Tags**

``{N}s``
  Exposure time tag (e.g., ``300s``, ``600s``). Used for matching darks to science frames.

``DARK_COEFF``
  Processed dark with fitted coefficients

``DARK_SYNTH``
  Synthetic dark from OBJECT frame


See Also
========

* :class:`~maroonx_instruments.maroonx.AstroDataMAROONX` - tag implementation
* `DRAGONS AstroData documentation <https://dragons.readthedocs.io/>`_ - general tag system