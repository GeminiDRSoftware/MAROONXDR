.. _tags:

****
Tags
****

Tags are metadata labels automatically assigned to MAROONX data files. They determine which recipes process the data and enable data selection for calibration matching.

MAROONX uses a **fiber-configuration-based tagging system** where tags are assigned based on the 5-fiber setup pattern rather than OBSTYPE keywords.


How MAROON-X Tags a File
=========================

When AstroData opens a MAROON-X file, it runs every
``@astro_data_tag``-decorated method on ``AstroDataMAROONX``, plus every
tag method inherited from ``AstroDataGemini`` and the AstroData core.
Each method inspects headers or extensions and either returns a tag set
or nothing. The file's final tag set is the union of what all these
methods return.

Some tag methods read the fiber illumination pattern, others read a
specific header keyword or the presence of a particular FITS extension.
The sections that follow document each category. The exact rule for
every tag lives in ``maroonx_instruments/maroonx/adclass.py``.


Fiber Configuration Patterns
=============================

The tag set for a MAROON-X file is determined primarily by the illumination pattern across its 5 fibers. The table below lists the 14 fiber configurations the pipeline recognises and the tags they produce. Additional tags (arm, exposure time, ``RAW`` / ``PROCESSED``, and instrument) are applied by separate rules and are documented below in `All Tags`_.

.. list-table::
   :header-rows: 1
   :widths: 15 45 40

   * - Fiber setup
     - Fiber pattern (1..5)
     - Tags produced
   * - ``FDDDF``
     - Flat lamp, Dark, Dark, Dark, Flat lamp
     - ``FLAT, CAL``
   * - ``DFFFD``
     - Dark, Flat lamp, Flat lamp, Flat lamp, Dark
     - ``FLAT, CAL``
   * - ``DDDDF``
     - Dark, Dark, Dark, Dark, Flat lamp
     - ``FLAT, CAL``
   * - ``DDDDE``
     - Dark, Dark, Dark, Dark, Etalon
     - ``DARK, CAL``
   * - ``SOOOE``
     - Sky, Target, Target, Target, Etalon
     - ``SCI, SPECT``
   * - ``DTTTT``
     - Dark, ThAr, ThAr, ThAr, ThAr
     - ``WAVECAL, SPECT, ThAr, CAL``
   * - ``DTTTE``
     - Dark, ThAr, ThAr, ThAr, Etalon
     - ``WAVECAL, SPECT, ThAr, CAL``
   * - ``DTTTD``
     - Dark, ThAr, ThAr, ThAr, Dark
     - ``WAVECAL, SPECT, ThAr, CAL``
   * - ``DEEEE``
     - Dark, Etalon, Etalon, Etalon, Etalon
     - ``WAVECAL, SPECT, ETALON, CAL``
   * - ``DEEEI``
     - Dark, Etalon, Etalon, Etalon, Iodine cell
     - ``WAVECAL, SPECT, ETALON, CAL``
   * - ``DLLLL``
     - Dark, LFC, LFC, LFC, LFC
     - ``WAVECAL, SPECT, LFC, CAL``
   * - ``DLLLE``
     - Dark, LFC, LFC, LFC, Etalon
     - ``WAVECAL, SPECT, LFC, CAL``
   * - ``DLLLD``
     - Dark, LFC, LFC, LFC, Dark
     - ``WAVECAL, SPECT, LFC, CAL``
   * - ``DEEEL``
     - Dark, Etalon, Etalon, Etalon, LFC
     - ``WAVECAL, SPECT, LFC, CAL``

.. note::

   **Additional tags applied per file.** Every file also carries ``MAROONX`` (instrument), one of ``BLUE`` / ``RED`` / ``BUNDLE`` (arm), and ``{N}s`` (exposure time). Raw files carry ``RAW``; processed files carry ``PROCESSED``. These are applied by separate tag methods (``_tag_instrument``, ``_tag_arm``, ``_tag_exptime``, ``_status_processed_maroonx_cals`` in ``adclass.py``), not by fiber configuration.

.. note::

   **Dark subtype tags.** A dark file (fiber setup ``DDDDE``) additionally carries ``DARK_COEFF`` if it has a ``COEFF_Z0`` extension (fitted dark coefficients), or ``DARK_SYNTH`` if its ``OBSTYPE`` is ``OBJECT`` (synthetic dark derived from a science exposure). These subtype tags are applied in the same ``_tag_dark`` method that assigns ``DARK, CAL``.




Tag to Recipe Mapping
======================

Each row below lists a tag set and the default recipe that ``reduce``
fires when it sees that tag set. These are the mappings a file must
satisfy to be processed automatically; any other recipe in the same
file must be requested explicitly with ``--recipe`` (or
``myreduce.recipename`` in the Python API).

.. list-table::
   :header-rows: 1
   :widths: 25 30 45

   * - Tags
     - Default recipe
     - Available recipes
   * - ``BUNDLE``
     - :meth:`~maroonxdr.maroonx.recipes.sq.recipes_BUNDLE.processBundle`
     - ``processBundle``
   * - ``DARK``
     - :meth:`~maroonxdr.maroonx.recipes.sq.recipes_DARK.makeProcessedDark`
     - ``makeProcessedDark``, ``makeDarkCoefficients``, ``testVARDark``,
       ``testRegressionDark``
   * - ``FLAT``
     - :meth:`~maroonxdr.maroonx.recipes.sq.recipes_FLAT_SPECT.makeProcessedFlat`
     - ``makeProcessedFlat``, ``makeProcessedFlatDFFFF``,
       ``makeStrayLightCheck``, ``makeFlatVarCheck``, ``measureBlaze``
   * - ``WAVECAL, ThAr``
     - :meth:`~maroonxdr.maroonx.recipes.sq.recipes_STATICWAVECAL.makeStaticWavecal`
     - ``makeStaticWavecal``
   * - ``WAVECAL``
     - :meth:`~maroonxdr.maroonx.recipes.sq.recipes_DYNAMIC_WAVECAL.makeDynamicWavecal`
     - ``makeDynamicWavecal``
   * - ``SCI``
     - :meth:`~maroonxdr.maroonx.recipes.sq.recipes_ECHELLE_SPECT.reduce`
     - ``reduce``, ``makeSyntheticDark``, ``applyBarycentricCorrection``,
       ``exportReducedBundle``, ``makeStripeExtractionCheck``

Non-default recipes are invoked explicitly with ``--recipe <name>`` on
the CLI or ``myreduce.recipename = '<name>'`` in the Python API. See the
:ref:`recipes` chapter for what each one does.



Using Tags for Data Selection
==============================

MAROON-X is an external DRAGONS instrument package (see :ref:`maroonxdr_prog_intro`),
so every ``typewalk`` and ``dataselect`` call must be told where to find
the AstroData tags with ``--adpkg maroonx_instruments``. Without that
flag, the MAROON-X tag methods never load and every file looks untagged.

``typewalk`` prints the tags of every FITS file in the current
directory:

.. code-block:: bash

    typewalk --adpkg maroonx_instruments -n

.. note::

   **Tags vs descriptors.** AstroData exposes header information in two
   ways. ``@astro_data_tag`` methods return labels for classification
   (``BLUE``, ``DARK``, ``PROCESSED``); ``@astro_data_descriptor`` methods
   return values (``ad.exposure_time()`` returns ``300.0``,
   ``ad.arm()`` returns ``'BLUE'``). ``dataselect --tags`` filters on
   the former, ``dataselect --expr`` on the latter, for example
   ``--expr 'exposure_time==300'``. The two can be combined.

``dataselect`` filters an explicit list of files by tag:

.. code-block:: bash

    # All wavelength calibrations
    dataselect --adpkg maroonx_instruments --tags WAVECAL *.fits

    # Blue arm science
    dataselect --adpkg maroonx_instruments --tags SCI,BLUE *.fits

    # Processed flats
    dataselect --adpkg maroonx_instruments --tags FLAT,PROCESSED *.fits

    # Raw darks with 300s exposure
    dataselect --adpkg maroonx_instruments --tags RAW,DARK,300s *.fits

Exclude tags with ``--xtags``:

.. code-block:: bash

    # Regular darks, excluding coefficients and synthetic
    dataselect --adpkg maroonx_instruments --tags DARK \
        --xtags DARK_COEFF,DARK_SYNTH *.fits

In Python scripts, importing ``maroonx_instruments`` registers the
MAROON-X tag methods; without that import ``select_data`` sees no
MAROON-X tags:

.. code-block:: python

    from gempy.adlibrary import dataselect
    import maroonx_instruments  # noqa: registers MAROON-X AstroData tags

    # Select blue arm science with 300s exposure
    sci_files = dataselect.select_data(
        all_files,
        tags=['BLUE', 'SCI', '300s']
    )


All Tags
=========

**Instrument/Observatory**

``GEMINI``, ``NORTH``, ``MAROONX``


**Arm Configuration** (mutually exclusive)

``BUNDLE`` : Raw GOA archive, blue and red arms as extensions

``BLUE`` : Blue arm data (450-670 nm)

``RED`` : Red arm data (650-920 nm)


**Data modality**

``SPECT`` : Spectroscopic data. Applied to every MAROON-X science and
calibration file.


**Purpose** (mutually exclusive)

``SCI`` : Science observation.

``CAL`` : Calibration frame. Applied to every non-science file.


**Calibration category** (mutually exclusive within ``CAL``)

``DARK`` : Dark frame.

``FLAT`` : Flat field.

``WAVECAL`` : Wavelength calibration. Further specialised by source:

  * ``ThAr`` : Thorium-argon arc lamp.
  * ``ETALON`` : Fabry-Perot etalon.
  * ``LFC`` : Laser frequency comb.

``BPM`` : Bad pixel mask. Applied when ``OBSTYPE=BPM``, not from fiber
setup.


**Processing Status**

``RAW`` : Unprocessed data

``PROCESSED`` : Processed by DRAGONS

``PREPARED`` : Basic preparation applied

**Special Tags**

``{N}s`` : Exposure time tag (e.g., ``300s``, ``600s``)

``DARK_COEFF`` : Processed dark with fitted coefficients

``DARK_SYNTH`` : Synthetic dark from science frame

``BARYCOR`` : Barycentric velocity correction recorded


See Also
========

* :class:`~maroonx_instruments.maroonx.AstroDataMAROONX` - tag implementation
* `DRAGONS AstroData documentation <https://dragons.readthedocs.io/>`_ - general tag system