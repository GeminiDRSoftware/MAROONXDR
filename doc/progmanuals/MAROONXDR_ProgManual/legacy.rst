.. _legacy:

*****************
Legacy References
*****************

This document provides a comprehensive comparison between the legacy MaroonX
HDF5-based pipeline products and the new DRAGONS FITS-based pipeline products.

File Structure Overview
=======================

The DRAGONS pipeline produces FITS files organized by extensions, while the
legacy pipeline used HDF5 format with hierarchical group structures. This
section maps legacy HDF5 datasets to their DRAGONS FITS extension equivalents.

Legacy File Naming Convention
------------------------------

Native MAROON-X raw files follow this naming pattern:

  ``YYYYMMDDTHHmmSSZ_SSSSS_C_nnnn.fits``. Example: ``20241115T123456Z_SOOOE_b_0600.fits``

Where:

  * ``YYYYMMDD`` - UT date of exposure start
  * ``HHmmSS`` - UT time of exposure start
  * ``SSSSS`` - Light source configuration for each of the five MAROON-X fibers:

    * ``D`` - Dark
    * ``F`` - Flat field
    * ``O`` - Object (stellar target)
    * ``S`` - Sky background (fiber 1 only)
    * ``E`` - Etalon comb
    * ``T`` - ThAr arc
    * ``I`` - Iodine cell
    * ``L`` - Laser Frequency Comb

  * ``C`` - Camera arm (``b`` = blue, ``r`` = red)
  * ``nnnn`` - Exposure time in integer seconds

**Common exposure types:**

  * ``SOOOE`` - Science exposure (Sky, Object, Etalon)
  * ``DDDDE`` - Dark frame (exposure > 30s) or drift tracking (< 30s)
  * ``DEEEE`` - Etalon wavelength calibration
  * ``DTTT?`` - ThAr arc calibration
  * ``DFFFD`` - Flat field (science fibers only)
  * ``FDDDF`` - Flat field (sky and calibration fibers)
  * ``DDDDF`` - Flat field (some FDDDF cases have fiber 1 un-lit)
  * ``DLLLL`` - Laser frequency comb calibration

Gemini Observatory Archive Naming convention
--------------------------------------------
  
  ``NYYYYMMDDMnnnn.fits``. Example: ``N20241115M0001.fits``

Where ``N`` = Gemini North, ``M`` = MAROON-X, and ``nnnn`` is a sequential counter
within each UT day. This is a MEF container for the BLUE (index 0) and RED (index 1) files, which can be restored
to the native format using the ``processBundle`` recipe. 

.. code-block:: pycon

   >>> import astrodata
   >>> import maroonx_instruments
   >>> ad = astrodata.open('N20241124M2838.fits')
   >>> ad.info()
   Filename: N20241124M2838.fits
   Tags: BUNDLE GEMINI MAROONX NORTH RAW SCI SPECT UNPREPARED

   Pixels Extensions
   Index  Content                  Type              Dimensions     Format
   [ 0]   science                  NDAstroData       (4400, 4400)   uint16
   [ 1]   science                  NDAstroData       (4400, 4400)   uint16

   Other Extensions
                  Type        Dimensions
   .EXPOSUREMETER Table       (278, 3)

.. note::
   The order of the arms in the bundle is ALWAYS guaranteed to be BLUE, ``ad[0]``, then RED, ``ad[1]``.


Dark Frame Products
===================

Raw Darks
---------

.. topic:: Legacy Format (FITS)
   :class: legacy-block

   Raw dark frames use single-extension FITS files, one per camera arm:

   * **Naming**: ``YYYYMMDDTHHmmSSZ_DDDDE_[b|r]_nnnn.fits``

     * ``DDDDE``: Fiber configuration (Dark-Dark-Dark-Dark-Etalon)
     * ``[b|r]``: Blue or red camera arm
     * ``nnnn``: Exposure time in seconds (e.g., ``0060``, ``1800``)

   * **Structure**: Single ImageHDU with 4400×4400 pixels (uint16)

.. topic:: DRAGONS Format (FITS/MEF)
   :class: dragons-block

   Raw dark frames from GOA (Gemini Observatory Archive) use multi-extension FITS to store both arms:

   * **Naming**: ``NYYYYMMDDMiiii.fits``

     * ``N``: Gemini North
     * ``M``: MAROON-X instrument code
     * ``iiii``: Sequential counter within UT day

   * **Structure**: Bundle format containing both arms:

      .. code-block:: pycon

         >>> import astrodata
         >>> import maroonx_instruments

         >>> ad = astrodata.open('N20241115M3455.fits')
         >>> ad.info()
         Filename: N20241115M3455.fits
         Tags: BUNDLE CAL DARK GEMINI MAROONX NORTH RAW SPECT UNPREPARED

         Pixels Extensions
         Index  Content                  Type              Dimensions     Format
         [ 0]   science                  NDAstroData       (4400, 4400)   uint16
         [ 1]   science                  NDAstroData       (4400, 4400)   uint16

         Other Extensions
                        Type        Dimensions
         .EXPOSUREMETER Table       (139, 3)

   * **Special Key Headers**:

     * ``INSTRUME='MAROON-X'``
     * ``OBSTYPE='DARK'``
     * ``FIBER[1-5]='Dark','Dark','Dark','Dark','Etalon'``
     * ``ORIGNAME``: Original GOA filename

   * **Processing**: Use ``processBundle()`` recipe to separate arms before reduction


Master Darks
------------

.. topic:: Legacy Format (FITS)
   :class: legacy-block

   * **Naming**: ``YYYYMMDDTHH_masterdark_mean_DDDDE_[b|r]_nnnn.fits``

   * **Structure**: Single ImageHDU with 4400×4400 pixels (float64)

   * **Standard Exposure Times**: 60, 120, 300, 600, 900, 1200, 1800 seconds
   * **Code Reference**: ``reduce/recipes/make_master_darks.py``

.. topic:: DRAGONS Format (FITS)
   :class: dragons-block

   * **Naming**: The name of the first input dark is taken, ``YYYYMMDDTHHmmSSZ_DDDDE_[b|r]_nnnn_dark.fits`` 
      (e.g., ``20241116T001751Z_DDDDE_r_1800_dark.fits``)
   
   * **Structure**: Single-arm file with 1 science extension:

      .. code-block:: pycon
         
         >>> import astrodata
         >>> import maroonx_instruments

         >>> ad = astrodata.open('20241116T001751Z_DDDDE_r_1800_dark.fits')
         >>> ad.info()
         Filename: 20241116T001751Z_DDDDE_r_1800_dark.fits
         Tags: 1800s CAL DARK GEMINI MAROONX NORTH OVERSCAN_SUBTRACTED PREPARED PROCESSED
            RED SPECT

         Pixels Extensions
         Index  Content                  Type              Dimensions     Format
         [ 0]   science                  NDAstroData       (4400, 4400)   float64
                  .variance             ADVarianceUncerta (4400, 4400)   float64
                  .mask                 ndarray           (4400, 4400)   uint16

         Other Extensions
                        Type        Dimensions
         .EXPOSUREMETER Table       (465, 3)
         .HISTORY       Table       (3, 4)
         .PROVENANCE    Table       (2, 4)


   * **Processing Recipe** (``recipes_DARK.py::makeProcessedDark``):

   * **Key Headers**:

      .. todo::
         List all relevant headers keywords created after reduction


Dark Coefficient Files
----------------------

.. topic:: Legacy Format (Numpy NPZ)
   :class: legacy-block

   Coefficient files enable synthetic dark generation as a function of exposure time:

   * **Naming**: ``masterdarks_coeffs_YYYYMMxx_[blue|red].npz``
   * **Format**: Numpy compressed archive (NPZ) with three arrays:

     * ``z0``: Slope coefficients (4400×4400, float64)
     * ``z1``: Intercept coefficients (4400×4400, float64)
     * ``logexptime``: log(exposure times) used for fitting

   * **Code Reference**: ``reduce/recipes/make_coeffs_from_masterdarks.py``

.. topic:: DRAGONS Format (FITS)
   :class: dragons-block

   Coefficient files store the same log-linear model using FITS table extensions:

   * **Naming**: ``YYYYMMDDTHHmmSSZ_DDDDE_[b|r]_nnnn_darkCoefficients.fits``
   * **Structure**: Single-arm with coefficient extensions

      .. code-block:: pycon

         >>> ad = astrodata.open('20241115T190028Z_DDDDE_b_0120_darkCoefficients.fits')
         >>> ad.info()
         Filename: 20241115T190028Z_DDDDE_b_0120_darkCoefficients.fits
         Tags: 120s BLUE CAL DARK DARK_COEFF GEMINI MAROONX NORTH OVERSCAN_SUBTRACTED
            PREPARED PROCESSED SPECT

         Pixels Extensions
         Index  Content                  Type              Dimensions     Format
         [ 0]   science                  NDAstroData       (1, 1)         float64
                  .variance             ADVarianceUncerta (4400, 4400)   float64
                  .mask                 ndarray           (4400, 4400)   uint16
                  .COEFF_Z0             ndarray           (4400, 4400)   float64
                  .COEFF_Z1             ndarray           (4400, 4400)   float64
                  .LOGEXPTIME           Table             (7, 3)         n/a

         Other Extensions
                        Type        Dimensions
         .EXPOSUREMETER Table       (139, 3)
         .HISTORY       Table       (4, 4)
         .PROVENANCE    Table       (4, 4)

   * **Processing Recipe** (``recipes_DARK.py::makeDarkCoefficients``):

   * **Tags**: ``{'GEMINI', '[BLUE|RED]', 'DARK', 'DARK_COEFF', 'CAL', 'PROCESSED'}``

Synthetic Darks
---------------

.. todo::
   Add description of the synthetic darks. Uses SOOOE frames

Flat Field Products
===================

Raw Flats
---------

.. topic:: Legacy Format (FITS)
   :class: legacy-block

   Raw flat field frames use single-extension FITS files, one per camera arm:

   * **Naming**: ``YYYYMMDDTHHmmSSZ_SSSSS_[b|r]_nnnn.fits``

     * ``SSSSS``: Fiber configuration (``DFFFD`` for science fibers 2,3,4; ``DDDDF`` for calibration fiber 5)
     * ``nnnn``: Exposure time in seconds

   * **Structure**: Single ImageHDU with 4400×4400 pixels (uint16)

.. topic:: DRAGONS Format (FITS/MEF)
   :class: dragons-block

   Raw flat fields from GOA use bundle format containing both arms:

   * **Naming**: ``NYYYYMMDDMiiii.fits``
   * **Structure**: Bundle format with EXPOSUREMETER extension.

      .. code-block:: pycon   

         >>>> ad = astrodata.open('N20241114M3300.fits')
         >>>> ad.info()
         Filename: N20241114M3300.fits
         Tags: BUNDLE CAL FLAT GEMINI MAROONX NORTH RAW SPECT UNPREPARED

         Pixels Extensions
         Index  Content                  Type              Dimensions     Format
         [ 0]   science                  NDAstroData       (4400, 4400)   uint16
         [ 1]   science                  NDAstroData       (4400, 4400)   uint16

         Other Extensions
                        Type        Dimensions
         .EXPOSUREMETER Table       (2239, 3)

   * **Processing**: Use ``processBundle()`` to separate arms before reduction


Master Flats
------------

.. topic:: Legacy Format (HDF5)
   :class: legacy-block

   The legacy pipeline produces a 4-stage flat workflow culminating in HDF5 format:

   * **Final Product**: ``YYYYMMDDTHH_masterflat_backgroundsubtracted_FFFFF_[b|r]_nnnn.hdf``
   * **Format**: HDF5 with nested structure:

     .. code-block:: text

        stripe_indices/
          fiber  (4072, 3954) int8                 # Pixel-to-fiber map
          order  (4072, 3954) int8                 # Pixel-to-order map
        box_extraction/fiber_N/order_M/            # 1D box-summed spectra
        optimal_extraction/fiber_N/order_M/        # Optimally extracted spectra
        optimal_var/fiber_N/order_M/               # Variance arrays
        wavelengths_static/fiber_N/order_M/        # Wavelength grids
        extraction_parameters/fiber_N/order_M/     # Extraction params

   * **Code References**: ``reduce/recipes/make_master_flats.py``, ``reduce/recipes/backgroundfit.py``, ``reduce/extraction.py``


.. topic:: DRAGONS Format (FITS)
   :class: dragons-block

   The DRAGONS pipeline takes the selected raw flats with ``DFFFD`` and ``DDDDF|FDDDF`` and produces a processed flat:

   * **Naming**: The names takes the name of the first raw flat in the input list and appends the suffix ``_[DFFFF|FFFFF]_flat``, e.g.,
      ``20241114T190714Z_DDDDF_b_0007_DFFFF_flat.fits``.
   * **Structure**: Each arm is stored individually.

      .. code-block:: pycon

         >>>> ad = astrodata.open('20241114T190714Z_DDDDF_b_0007_DFFFF_flat.fits')
         >>>> ad.info()
         Filename: 20241114T190714Z_DDDDF_b_0007_DFFFF_flat.fits
         Tags: 7s BLUE CAL FLAT GEMINI MAROONX NORTH OVERSCAN_SUBTRACTED OVERSCAN_TRIMMED
            PREPARED PROCESSED SPECT

         Pixels Extensions
         Index  Content                  Type              Dimensions     Format
         [ 0]   science                  NDAstroData       (4072, 3954)   float64
                  .variance             ADVarianceUncerta (4072, 3954)   float64
                  .mask                 ndarray           (4072, 3954)   uint16
                  .INDEX_FIBER          ndarray           (4072, 3954)   int64
                  .INDEX_ORDER          ndarray           (4072, 3954)   int64
                  .REMOVED_STRIPES      ndarray           (4, 6)         float64
                  .STRIPES_FIBERS       ndarray           (4,)           int64
                  .STRIPES_ID           Table             (24, 34)       n/a

         Other Extensions
                        Type        Dimensions
         .EXPOSUREMETER Table       (118, 3)
         .HISTORY       Table       (4, 4)
         .PROVENANCE    Table       (2, 4)

   * **Processing Recipe**: ``recipes_FLAT_SPECT.py::makeProcessedFlat``


Wavelength Calibration Products
===============================

Raw Wavelength Calibrations
---------------------------

.. topic:: Legacy Format (FITS)
   :class: legacy-block

   Raw wavelength calibration frames use single-extension FITS files, one per camera arm:

   * **Naming**: ``YYYYMMDDTHHmmSSZ_SSSSS_[b|r]_nnnn.fits``
      Where the fiber setup ``SSSSS`` is ``DEEEE`` (etalon), or ``DLLLL?`` (LFC)

   * **Structure**: Single ImageHDU with 4400×4400 pixels (uint16)

.. topic:: DRAGONS Format (FITS/MEF)
   :class: dragons-block

   Raw wavelength calibrations from GOA use bundle format containing both arms:

   * **Naming**: ``NYYYYMMDDMiiii.fits``
   * **Structure**: Bundle format with EXPOSUREMETER extension

      .. code-block:: pycon

         >>>> ad = astrodata.open('N20241124M0554.fits')
         >>>> ad.info()
         Filename: N20241124M0554.fits
         Tags: BUNDLE CAL GEMINI LFC MAROONX NORTH RAW SPECT UNPREPARED WAVECAL

         Pixels Extensions
         Index  Content                  Type              Dimensions     Format
         [ 0]   science                  NDAstroData       (4400, 4400)   uint16
         [ 1]   science                  NDAstroData       (4400, 4400)   uint16

         Other Extensions
                        Type        Dimensions
         .EXPOSUREMETER Table       (2057, 3)

   * **Tags**: Wavelength calibration type identified by tag combinations:

     * ``{'WAVECAL', 'ETALON'}`` - Etalon calibration frames (DEEEE pattern)
     * ``{'WAVECAL', 'LFC'}`` - Laser frequency comb (DLLLL? pattern)
     * ``{'WAVECAL', 'ThAr'}`` - Thorium-argon arc lamp (DTTT? pattern)

   * **Processing**: Use ``processBundle()`` to separate arms before wavelength calibration


Static Wavelength Solutions
---------------------------

.. topic:: Legacy Format (HDF5)
   :class: legacy-block

   Static wavelength solutions provide initial calibration:

   * **Naming**: ``20241124T030227Z_DEEEE_b_0030.hdf``
   * **Format**: HDF5 with nested structure:

       .. code-block:: text
   
         /wavelengths_static/
            fiber_1/             # Empty
            fiber_2/ ...         # M orders
            fiber_3/ ...         # M orders
            fiber_4/ ...         # M orders
            fiber_5/ ...         # M orders


.. topic:: DRAGONS Format (FITS)
   :class: dragons-block
  
   This Static Wavelength solution is loaded to a wavecal file in the ``staticWavelengthSolution()`` primitive.
   
   * **Storage in Processed Files**: Stored as ``WLS_STATIC_FIBER_*`` extensions

      .. code-block:: pycon

         >>>> ad = astrodata.open('20241124T030227Z_DEEEE_b_0030_wavecal.fits')
         >>>> ad.info()
         Filename: 20241124T030227Z_DEEEE_b_0030_wavecal.fits
         Tags: 30s BLUE CAL ETALON GEMINI MAROONX NORTH OVERSCAN_SUBTRACTED
            OVERSCAN_TRIMMED PREPARED PROCESSED SPECT WAVECAL

         Pixels Extensions
         Index  Content                  Type              Dimensions     Format
         [ 0]   science                  NDAstroData       (4072, 3954)   float32
                  ...
                  .WLS_STATIC_FIBER_1   ndarray           (1, 1)         float64
                  .WLS_STATIC_FIBER_2   ndarray           (34, 3954)     float64
                  .WLS_STATIC_FIBER_3   ndarray           (34, 3954)     float64
                  .WLS_STATIC_FIBER_4   ndarray           (34, 3954)     float64
                  .WLS_STATIC_FIBER_5   ndarray           (34, 3954)     float64


Dynamic Wavelength Solutions
-----------------------------

.. topic:: Legacy Format (HDF5)
   :class: legacy-block

   Dynamic etalon-based solutions:

   * **Processing Script**: ``analyze/recipes/batch_etalon_spline_wls.py``
   * **File structure**: In addition to ``wavelengths_static``, the HDF5 file contains

     .. code-block:: text

        /wavelengths_dynamic/
          fiber_2/
            M orders: float64    # M wavelength arrays [nm]
          fiber_3/ ...
          fiber_4/ ...
          fiber_5/ ...

        /etalon_peak_parameters/
          peaks/           # PyTables Table with Gaussian fit parameters
          polynomials/     # Polynomial coefficients

        /peak_data/        # Pandas DataFrame
          # Multi-index: (fiber, order, m)
          # Columns: center, amplitude, sigma, m, m_fraction,
          #          wavelength_by_thar, dispersion_mps

        /wavelengths/      # Soft link → /wavelengths_dynamic/

   .. todo::
      Add instrument drift keywords here

.. topic:: DRAGONS Format (FITS)
   :class: dragons-block

   Dynamic wavelength solutions stored as 2D arrays in processed wavecal files:

   * **Naming**: ``YYYYMMDDTHHmmSSZ_DEEEE_[b|r]_nnnn_wavecal.fits``
   * **Recipe**: ``makeDynamicWavecal`` (``recipes_DYNAMIC_WAVECAL.py``)
   * **Structure**: Single-arm file with wavelength solution extensions:

      .. code-block:: pycon

         >>> ad = astrodata.open('20241124T030227Z_DEEEE_r_0004_wavecal.fits')
         >>> ad.info()
         Filename: 20241124T030227Z_DEEEE_r_0004_wavecal.fits
         Tags: RED ETALON 4s NORTH PREPARED OVERSCAN_TRIMMED WAVECAL MAROONX
               OVERSCAN_SUBTRACTED GEMINI SPECT PROCESSED CAL

         Pixels Extensions
         Index  Content                  Type              Dimensions     Format
         [ 0]   science                  NDAstroData       (4080, 4036)   float32
                  .WLS_STATIC_FIBER_1   ndarray           (1, 1)         float64
                  .WLS_STATIC_FIBER_2   ndarray           (28, 4036)     float64
                  .WLS_STATIC_FIBER_3   ndarray           (28, 4036)     float64
                  .WLS_STATIC_FIBER_4   ndarray           (28, 4036)     float64
                  .WLS_STATIC_FIBER_5   ndarray           (28, 4036)     float64
                  .WLS_DYNAMIC_FIBER_1  ndarray           (1, 1)         float64
                  .WLS_DYNAMIC_FIBER_2  ndarray           (28, 4036)     float64
                  .WLS_DYNAMIC_FIBER_3  ndarray           (28, 4036)     float64
                  .WLS_DYNAMIC_FIBER_4  ndarray           (28, 4036)     float64
                  .WLS_DYNAMIC_FIBER_5  ndarray           (28, 4036)     float64
                  .REDUCED_ORDERS_FIBER_2 ndarray         (28,)          float64
                  .PEAKS                Table             (43710, 11)    n/a
                  .PEAK_DATA            Table             (43646, 16)    n/a
                  .POLY                 Table             (112, 10)      n/a
                  ... (other extensions)


   * **Etalon Fitting Results**:

     * ``PEAKS``: Table with Gaussian fit parameters (11 columns)
     * ``PEAK_DATA``: Refined peaks with wavelength assignments (16 columns)
     * ``POLY``: Polynomial parameters for peak profile modeling (10 columns)

   * **Processing Primitives**:

     1. ``getPeaksAndPolynomials()`` - Iterative etalon peak fitting
     2. ``staticWavelengthSolution()`` - Load static calibration
     3. ``fitAndApplyEtalonWls()`` - Fit dynamic solution with cubic splines


Simultaneous Wavelength Solutions
----------------------------------

.. topic:: Legacy Format (HDF5)
   :class: legacy-block

   Simultaneous etalon drift correction applied to science frames:

   * **Processing Script**: ``analyze/recipes/batch_science_spline_wls_dynamic.py``
   * **Target Files**: SOOOE science frames with simultaneous etalon (fiber 5)

   * **HDF5 Structure**:

     .. code-block:: text

        /wavelengths_simultaneous/
          fiber_2/
            M orders: float64    # Science fiber with drift correction
            ...
          fiber_3/ ...
          fiber_4/ ...
          fiber_5/               # Reference etalon fiber
            M orders: float64
            ...

        /wavelengths/            # Soft link → /wavelengths_simultaneous/


.. topic:: DRAGONS Format (FITS)
   :class: dragons-block

   Simultaneous wavelength solutions applied via science reduction recipe:

   * **Naming**: ``YYYYMMDDTHHmmSSZ_SOOOE_[b|r]_nnnn_reduced.fits``
   * **Recipe**: ``reduce`` in ``recipes_ECHELLE_SPECT.py``
   * **Structure**: Science frame with simultaneous wavelength extensions:

      .. code-block:: pycon

         >>> ad = astrodata.open('20241124T041907Z_SOOOE_r_0300_reduced.fits')
         >>> ad.info()
         Filename: 20241124T041907Z_SOOOE_r_0300_reduced.fits
         Tags: 300s GEMINI MAROONX NORTH OVERSCAN_SUBTRACTED OVERSCAN_TRIMMED PREPARED
            PROCESSED PROCESSED_SCIENCE RED SCI SPECT

         Pixels Extensions
         Index  Content                  Type              Dimensions     Format
         [ 0]   science                  NDAstroData       (4080, 4036)   float32
                  .WLS_SIMULTANEOUS_FIB ndarray           (1, 1)         float64
                  .WLS_SIMULTANEOUS_FIB ndarray           (28, 4036)     float64
                  .WLS_SIMULTANEOUS_FIB ndarray           (28, 4036)     float64
                  .WLS_SIMULTANEOUS_FIB ndarray           (28, 4036)     float64
                  .WLS_SIMULTANEOUS_FIB ndarray           (28, 4036)     float64

     The extensions name appear cut, they correspond to each fiber:
     * ``WLS_SIMULTANEOUS_FIBER_1`` - Empty (sky fiber)
     * ``WLS_SIMULTANEOUS_FIBER_2`` - Drift-corrected wavelength for science fiber 2
     * ``WLS_SIMULTANEOUS_FIBER_3`` - Drift-corrected wavelength for science fiber 3
     * ``WLS_SIMULTANEOUS_FIBER_4`` - Drift-corrected wavelength for science fiber 4
     * ``WLS_SIMULTANEOUS_FIBER_5`` - Reference etalon wavelength

   * **Key Primitive**: ``applyWavelengthSolution()``

Science Frame Products
======================

Raw Science Frames
------------------

.. topic:: Legacy Format (FITS)
   :class: legacy-block

   Raw science frames use single-extension FITS files, one per camera arm:

   * **Naming**: ``YYYYMMDDTHHmmSSZ_SOOOE_[b|r]_nnnn.fits``

     * ``SOOOE``: Fiber configuration (Sky-Object-Object-Object-Etalon)
     * ``[b|r]``: Blue or red camera arm
     * ``nnnn``: Exposure time in seconds

   * **Structure**: Single ImageHDU with 4400×4400 pixels (uint16)

.. topic:: DRAGONS Format (FITS/MEF)
   :class: dragons-block

   Raw science frames from GOA use bundle format containing both arms:

   * **Naming**: ``NYYYYMMDDMiiii.fits``
   * **Structure**: Bundle format with EXPOSUREMETER extension

      .. code-block:: pycon

         >>> ad = astrodata.open('N20241124M2838.fits')
         >>> ad.info()
         Filename: N20241124M2838.fits
         Tags: BUNDLE GEMINI MAROONX NORTH RAW SCI SPECT UNPREPARED

         Pixels Extensions
         Index  Content                  Type              Dimensions     Format
         [ 0]   science                  NDAstroData       (4400, 4400)   uint16
         [ 1]   science                  NDAstroData       (4400, 4400)   uint16

         Other Extensions
                        Type        Dimensions
         .EXPOSUREMETER Table       (278, 3)
         >>>> ad.fiber_setup()
         ['Sky', 'Target', 'Target', 'Target', 'Etalon']

   * **Processing**: Use ``processBundle()`` to separate arms before reduction


Reduced Science Frames
-----------------------

.. topic:: Legacy Format (HDF5/Pandas)
   :class: legacy-block

   The legacy pipeline produces final reduced science frames as HDF5/Pandas stores:

   * **Naming**: ``YYYYMMDDTHHmmSSZ_SOOOE_x_nnnn.hd5``

     * ``x``: Indicates combined arm data (not ``b`` or ``r``)
     * ``nnnn``: Exposure time in seconds

   * **Format**: HDF5 file containing Pandas DataFrames with multi-index structure:

     .. code-block:: text

        /spec_blue/           # Pandas DataFrame
          Index: [Fiber, Order]
          Columns: ['box_extraction', 'wavelengths', 'optimal_extraction', 'optimal_var']

        /spec_red/            # Pandas DataFrame
          Index: [Fiber, Order]
          Columns: ['box_extraction', 'wavelengths', 'optimal_extraction', 'optimal_var']

        /blaze_blue/          # Blaze function DataFrame
        /blaze_red/           # Blaze function DataFrame
        /header_blue/         # FITS header as DataFrame
        /header_red/          # FITS header as DataFrame

   * **Code References**:

     * ``analyze/recipes/combine_science_fibers.py`` - Fiber combination that adds Fiber 6|7
     * ``analyze/maroonxspectrum.py`` - HDF5 spectrum class


.. topic:: DRAGONS Format (FITS)
   :class: dragons-block

   DRAGONS pipeline produces reduced science frames as multi-extension FITS files:

   * **Naming**: ``YYYYMMDDTHHmmSSZ_SOOOE_[b|r]_nnnn_reduced.fits``

     * Separate files for blue and red arms
     * ``_reduced`` suffix indicates fully reduced science frame

   * **Structure**: Single-arm file with extracted spectra and wavelength solutions:

      .. code-block:: pycon

         >>> ad = astrodata.open('20241124T041907Z_SOOOE_r_0300_reduced.fits')
         >>> ad.info()
         Filename: 20241124T041907Z_SOOOE_r_0300_reduced.fits
         Tags: 300s GEMINI MAROONX NORTH OVERSCAN_SUBTRACTED OVERSCAN_TRIMMED PREPARED
            PROCESSED PROCESSED_SCIENCE RED SCI SPECT

         Pixels Extensions
         Index  Content                  Type              Dimensions     Format
         [ 0]   science                  NDAstroData       (4080, 4036)   float32
                  .variance             ADVarianceUncerta (4080, 4036)   float64
                  .mask                 ndarray           (4080, 4036)   uint16

                  # Box Extraction (5 fibers)
                  .BOX_REDUCED_FIBER_1  ndarray           (1, 1)         float64
                  .BOX_REDUCED_FIBER_2  ndarray           (28, 4036)     float64
                  .BOX_REDUCED_FIBER_3  ndarray           (28, 4036)     float64
                  .BOX_REDUCED_FIBER_4  ndarray           (28, 4036)     float64
                  .BOX_REDUCED_FIBER_5  ndarray           (28, 4036)     float64

                  # Box Extraction Variance
                  .BOX_REDUCED_VAR_1    ndarray           (1, 1)         float64
                  .BOX_REDUCED_VAR_2    ndarray           (28, 4036)     float64
                  .BOX_REDUCED_VAR_3    ndarray           (28, 4036)     float64
                  .BOX_REDUCED_VAR_4    ndarray           (28, 4036)     float64
                  .BOX_REDUCED_VAR_5    ndarray           (28, 4036)     float64

                  # Flat-fielded Spectra
                  .BOX_REDUCED_FLAT_1   ndarray           (1, 1)         float64
                  .BOX_REDUCED_FLAT_2   ndarray           (28, 4036)     float64
                  .BOX_REDUCED_FLAT_3   ndarray           (28, 4036)     float64
                  .BOX_REDUCED_FLAT_4   ndarray           (28, 4036)     float64
                  .BOX_REDUCED_FLAT_5   ndarray           (28, 4036)     float64

                  # Optimal Extraction (5 fibers + combined fiber 6)
                  .OPTIMAL_REDUCED_FIBER_1  ndarray       (1, 1)         float64
                  .OPTIMAL_REDUCED_FIBER_2  ndarray       (28, 4036)     float64
                  .OPTIMAL_REDUCED_FIBER_3  ndarray       (28, 4036)     float64
                  .OPTIMAL_REDUCED_FIBER_4  ndarray       (28, 4036)     float64
                  .OPTIMAL_REDUCED_FIBER_5  ndarray       (1, 1)         float64
                  .OPTIMAL_REDUCED_FIBER_6  ndarray       (28, 4036)     float64

                  # Optimal Extraction Variance
                  .OPTIMAL_REDUCED_VAR_1    ndarray       (1, 1)         float64
                  .OPTIMAL_REDUCED_VAR_2    ndarray       (28, 4036)     float64
                  .OPTIMAL_REDUCED_VAR_3    ndarray       (28, 4036)     float64
                  .OPTIMAL_REDUCED_VAR_4    ndarray       (28, 4036)     float64
                  .OPTIMAL_REDUCED_VAR_5    ndarray       (1, 1)         float64
                  .OPTIMAL_REDUCED_VAR_6    ndarray       (28, 4036)     float64

                  # Wavelength Solutions (Simultaneous)
                  .WLS_SIMULTANEOUS_FIBER_1 ndarray       (1, 1)         float64
                  .WLS_SIMULTANEOUS_FIBER_2 ndarray       (28, 4036)     float64
                  .WLS_SIMULTANEOUS_FIBER_3 ndarray       (28, 4036)     float64
                  .WLS_SIMULTANEOUS_FIBER_4 ndarray       (28, 4036)     float64
                  .WLS_SIMULTANEOUS_FIBER_5 ndarray       (28, 4036)     float64
                  .WLS_SIMULTANEOUS_FIBER_6 ndarray       (28, 4036)     float64

                  # Static Wavelength Solutions
                  .WLS_STATIC_FIBER_1   ndarray           (1, 1)         float64
                  .WLS_STATIC_FIBER_2   ndarray           (28, 4036)     float64
                  .WLS_STATIC_FIBER_3   ndarray           (28, 4036)     float64
                  .WLS_STATIC_FIBER_4   ndarray           (28, 4036)     float64
                  .WLS_STATIC_FIBER_5   ndarray           (28, 4036)     float64

                  # Bad Pixel Masks per Fiber
                  .BPM_FIBER_1          ndarray           (1, 1)         float64
                  .BPM_FIBER_2          ndarray           (28, 4036)     int64
                  .BPM_FIBER_3          ndarray           (28, 4036)     int64
                  .BPM_FIBER_4          ndarray           (28, 4036)     int64
                  .BPM_FIBER_5          ndarray           (28, 4036)     int64

                  # Order Information
                  .REDUCED_ORDERS_FIBER_1   ndarray       (1, 1)         float64
                  .REDUCED_ORDERS_FIBER_2   ndarray       (28,)          float64
                  .REDUCED_ORDERS_FIBER_3   ndarray       (28,)          float64
                  .REDUCED_ORDERS_FIBER_4   ndarray       (28,)          float64
                  .REDUCED_ORDERS_FIBER_5   ndarray       (28,)          float64
                  .REDUCED_ORDERS_FIBER_6   ndarray       (28,)          float64

                  # Etalon Fitting Results (from simultaneous calibration)
                  .PEAKS                Table             (10887, 11)    n/a
                  .POLY                 Table             (28, 10)       n/a
                  .STRIPES_ID           Table             (24, 28)       n/a

         Other Extensions
                        Type        Dimensions
         .EXPOSUREMETER Table       (389, 3)
         .HISTORY       Table       (4, 4)
         .PROVENANCE    Table       (3, 4)

   * **Processing Recipe**: ``reduce`` in ``recipes_ECHELLE_SPECT.py``

   * **Combined Fiber**: Fiber 6 is the combination of fibers 2, 3, 4 (science fibers)


Extraction Products Comparison
-------------------------------

The following table maps legacy HDF5 data products to DRAGONS FITS extensions:

.. list-table:: Science Frame Data Product Mapping
   :header-rows: 1
   :widths: 35 35 30

   * - Legacy HDF5 Path
     - DRAGONS FITS Extension
     - Description
   * - ``/spec_[blue|red]/box_extraction``
     - ``.BOX_REDUCED_FIBER_N``
     - Box-summed 1D spectra
   * - ``/spec_[blue|red]/optimal_extraction``
     - ``.OPTIMAL_REDUCED_FIBER_N``
     - Optimally extracted spectra
   * - ``/spec_[blue|red]/optimal_var``
     - ``.OPTIMAL_REDUCED_VAR_N``
     - Variance for optimal extraction
   * - ``/spec_[blue|red]/wavelengths``
     - ``.WLS_SIMULTANEOUS_FIBER_N``
     - Wavelength calibration [nm]
   * - N/A (computed inline)
     - ``.BOX_REDUCED_VAR_N``
     - Box extraction errors
   * - N/A (computed inline)
     - ``.BOX_REDUCED_FLAT_N``
     - Flat-fielded box spectra
   * - ``fiber_6`` (combined)
     - ``.OPTIMAL_REDUCED_FIBER_6``
     - Combined fibers 2+3+4

.. note::
   * **Array Dimensions**:

     * Legacy: 1D arrays per order accessed via ``df.loc[fiber, order]``
     * DRAGONS: 2D arrays ``(N_orders, N_pixels)`` in single extension

   * **Fiber Indexing**: Both use fibers 2-4 for science, fiber 6 for combination
   * **Sky Fiber**: Fiber 1 has placeholder extensions (1×1 arrays) in DRAGONS
   

Configuration Files
===================

Static Wavelength Solutions
----------------------------

.. topic:: Legacy Static Wavelength Model (HDF5)
   :class: legacy-block

   Static wavelength solutions provide initial calibration accurate to ~500 m/s:

   * **Reference Parameter File**: ``wl_combined_final_etalon_peakmodel_2020.hdf``
   * **Location**: ``MaroonX_spectra_reduced/Maroonx_wls/202005xx/``
   * **Format**: HDF5 with nested structure:

     .. code-block:: text

        /wls_red/
          fiber_N/
            maxx: 4036
            x_norm: float64        # Normalized pixel coordinates
            orders: int            # Order numbers
            weights: float64       # Line weights
            wavelengths: float64   # Wavelength values [nm]
            poly_deg_x: 5
            poly_deg_y: 5

        /wls_blue/
          fiber_N/
            maxx: 3954
            x_norm: float64
            orders: float64
            weights: float64
            wavelengths: float64
            poly_deg_x: 5
            poly_deg_y: 5

        /dispersion/
          parameter: JSON-encoded lmfit Parameters
            # Etalon physical parameters:
            l, theta, n, disp_0...N
            

.. topic:: DRAGONS Static Wavelength Model (FITS)
   :class: dragons-block

   Static wavelength solutions stored as pre-computed lookup tables:

   * **Naming**: ``WLSTAT_[b|r].fits`` (blue/red arm)
   * **Location**: ``maroonxdr/maroonx/lookups/WLS/``
   * **Structure**: FITS file with extensions per fiber:

      .. code-block:: pycon

         >>>> ad = astrodata.open('maroonxdr/maroonx/lookups/WLS/WLSTAT_b.fits')
         >>>> ad.info()
         Filename: WLSTAT_b.fits
         Tags:

         Other Extensions
                        Type        Dimensions
         .FIBER_1       Table       (3954, 34)
         .FIBER_2       Table       (3954, 34)
         .FIBER_3       Table       (3954, 34)
         .FIBER_4       Table       (3954, 34)
         .FIBER_5       Table       (3954, 34)

      Each fiber extension is an astropy Table with wavelength values per order.

   * **Code reference**: ``utils/ref_wls_hdf2fits.py``


Bad Pixel Masks
---------------

.. topic:: Legacy Format (HDF5)
   :class: legacy-block

   Bad pixel masks are embedded in per-night configuration files:

   * **Location**: ``{DATAX}/MaroonX_spectra_reduced/Maroonx_configfiles/202411xx/config_{arm}.hdf``
   * **File Names**: ``config_[b|r].hdf``
   * **Format**: HDF5 files (~13.4 MB per file)

   * **HDF5 Structure**:

     .. code-block:: text

        /bad_pixel_map              (4400×4400, int64) - bad pixel locations
        /valid/                     (group) - valid pixel regions
        /wavelengths/               (group) - per-order wavelength grids
        /wavelengths_static/        (group) - initial wavelength solution
        /find_stripes/              (group) - stripe finding parameters
        /identify_stripes/          (group) - stripe identification parameters
        /correct_image_orientation/ (group) - orientation correction settings
        /used_thar_lines/           (group) - list of echelle orders with ThAr lines

.. topic:: DRAGONS Format (FITS)
   :class: dragons-block

   Bad pixel masks stored as standalone FITS files:

   * **Location**: ``maroonxdr/maroonx/lookups/BPM/``
   * **File Names**: ``BPM_[b|r]_0000.fits``
   * **Format**: FITS ImageHDU

   * **Structure**:

     .. code-block:: pycon

        >>> from astropy.io import fits
        >>> hdul = fits.open('maroonxdr/maroonx/lookups/BPM/BPM_b_0000.fits')
        >>> hdul.info()
        Filename: maroonxdr/maroonx/lookups/BPM/BPM_b_0000.fits
        No.    Name      Ver    Type      Cards   Dimensions   Format
        0  PRIMARY       1 PrimaryHDU       8   (4400, 4400)   int64



Stripe Identification Tables
-----------------------------

.. topic:: Legacy Format (HDF5)
   :class: legacy-block

   Stripe positions reference are stored in config HDF5 files:

   * **Location**: Within ``config_[b|r].hdf`` files
   * **HDF5 Groups**:

     * ``/find_stripes/`` - Stripe finding parameters
     * ``/identify_stripes/`` - Stripe identification parameters


.. topic:: DRAGONS Format (FITS)
   :class: dragons-block

   Stripe identification stored as standalone FITS binary tables:

   * **Location**: ``maroonxdr/maroonx/lookups/SID/``
   * **File Names**: ``SID_[b|r].fits``
   * **Format**: FITS BinTableHDU

   * **Structure**:

     .. code-block:: pycon

        >>> from astropy.io import fits
        >>> hdul = fits.open('maroonxdr/maroonx/lookups/SID/SID_b.fits')
        >>> hdul.info()
        Filename: SID_b.fits
        No.    Name      Ver    Type      Cards   Dimensions   Format
          0  PRIMARY       1 PrimaryHDU       8   ()
          1  SID           1 BinTableHDU     20   170R x 3C   [K, K, K]

     * Columns:
       * ``identify_fiber`` (int): Fiber number (1-5)
       * ``fiber_order`` (int): Echelle order number
       * ``fiber_position`` (int): Y-position on detector (pixels)





Cheat Sheet
===========

.. todo::
   Provide a quick reference table mapping legacy products to DRAGONS equivalents.
   Focus on array names within the files and key header keywords.