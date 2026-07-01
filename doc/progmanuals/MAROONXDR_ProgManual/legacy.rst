.. _legacy:

*****************
Legacy References
*****************

This document provides a comprehensive comparison between the legacy MaroonX
HDF5-based pipeline products and the new DRAGONS FITS-based pipeline products.

.. contents:: Contents
   :local:
   :depth: 1

File Structure Overview
=======================

The DRAGONS pipeline produces FITS files organized by extensions, while the
legacy pipeline used HDF5 format with hierarchical group structures. This
section maps legacy HDF5 datasets to their DRAGONS FITS extension equivalents.

Legacy File Naming Convention
------------------------------

Native MAROON-X raw files follow this naming pattern:

  ``YYYYMMDDTHHmmSSZ_SSSSS_C_nnnn.fits``. Example: ``20250717T144308Z_SOOOE_b_0300.fits``

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
  
  ``NYYYYMMDDMnnnn.fits``. Example: ``N20250717M5299.fits``

Where ``N`` = Gemini North, ``M`` = MAROON-X, and ``nnnn`` is a sequential counter
within each UT day. This is a MEF container for the BLUE (index 0) and RED (index 1) files, which can be restored
to the native format using the ``processBundle`` recipe. 

.. code-block:: pycon

   >>> import astrodata
   >>> import maroonx_instruments
   >>> ad = astrodata.open('N20250717M5299.fits')
   >>> ad.info()
   Filename: N20250717M5299.fits
   Tags: 300s BUNDLE GEMINI MAROONX NORTH RAW SCI SPECT UNPREPARED

   Pixels Extensions
   Index  Content                  Type              Dimensions     Format
   [ 0]   science                  NDAstroData       (4400, 4400)   uint16
   [ 1]   science                  NDAstroData       (4400, 4400)   uint16

   Other Extensions
                  Type        Dimensions
   .EXPOSUREMETER Table       (1027, 3)

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
     * ``nnnn``: Exposure time in integer seconds

   * **Example**: ``20250817T182823Z_DDDDE_b_0600.fits``
   * **Structure**: Single-extension ``PrimaryHDU`` with a ``(4400, 4400)`` ``uint16`` array.

.. topic:: DRAGONS Format (FITS/MEF bundle)
   :class: dragons-block

   Raw dark frames from GOA (Gemini Observatory Archive) use a multi-extension FITS bundle that stores both arms:

   * **Naming**: ``NYYYYMMDDMiiii.fits``

     * ``N``: Gemini North
     * ``M``: MAROON-X instrument code
     * ``iiii``: Sequential counter within UT day

   * **Example**: ``N20250707M6246.fits`` (300 s DARK, ``DDDDE`` fiber configuration)
   * **Structure**: Bundle format containing both arms:

      .. code-block:: pycon

         >>> import astrodata
         >>> import maroonx_instruments
         >>> ad = astrodata.open('N20250707M6246.fits')
         >>> ad.info()
         Filename: N20250707M6246.fits
         Tags: 300s BUNDLE CAL DARK GEMINI MAROONX NORTH RAW SPECT UNPREPARED

         Pixels Extensions
         Index  Content                  Type              Dimensions     Format
         [ 0]   science                  NDAstroData       (4400, 4400)   uint16
         [ 1]   science                  NDAstroData       (4400, 4400)   uint16

         Other Extensions
                        Type        Dimensions
         .EXPOSUREMETER Table       (69, 3)

   * **Tags**: ``['300s', 'BUNDLE', 'CAL', 'DARK', 'GEMINI', 'MAROONX', 'NORTH', 'RAW', 'SPECT', 'UNPREPARED']``
   * **Processing**: Use the ``processBundle`` recipe to separate the two arms into single-arm files before reduction.


Master Darks
------------

.. topic:: Legacy Format (FITS)
   :class: legacy-block

   * **Example**: ``20250707T17_masterdark_mean_DDDDE_b_0300.fits`` (300 s blue)
   * **Naming**: ``YYYYMMDDTHH_masterdark_mean_DDDDE_[b|r]_nnnn.fits``
   * **Structure**: single-extension ``PrimaryHDU`` holding a ``(4400, 4400)`` float64 array.
   * **Standard exposure times**: 60, 120, 300, 600, 900, 1200, 1800 seconds.

   * **Code reference**: ``reduce/recipes/make_master_darks.py``

.. topic:: DRAGONS Format (FITS)
   :class: dragons-block

   * **Example**: ``20250707T172105Z_DDDDE_b_0300_dark.fits`` (300 s blue)
   * **Naming**: takes the filename of the first input dark and appends ``_dark``.

   * **Structure**: Master darks are ``OVERSCAN_SUBTRACTED`` but not ``OVERSCAN_TRIMMED``,
     so the science array keeps the full ``(4400, 4400)`` shape. Master flats, wavecals, and
     reduced science files are trimmed to ``(4072, 3954)``; darks are not.

      .. code-block:: pycon

         >>> import astrodata
         >>> import maroonx_instruments

         >>> ad = astrodata.open('20250707T172105Z_DDDDE_b_0300_dark.fits')
         >>> ad.info()
         Filename: 20250707T172105Z_DDDDE_b_0300_dark.fits
         Tags: 300s BLUE CAL DARK GEMINI MAROONX NORTH OVERSCAN_SUBTRACTED PREPARED
             PROCESSED SPECT

         Pixels Extensions
         Index  Content                  Type              Dimensions     Format
         [ 0]   science                  NDAstroData       (4400, 4400)   float32
                   .variance             ADVarianceUncerta (4400, 4400)   float32
                   .mask                 ndarray           (4400, 4400)   uint16

         Other Extensions
                        Type        Dimensions
         .EXPOSUREMETER Table       (69, 3)
         .HISTORY       Table       (2, 4)
         .PROVENANCE    Table       (2, 4)

   * **Processing recipe**: ``recipes_DARK.py::makeProcessedDark``


Dark Coefficient Files
----------------------

.. topic:: Legacy Format (Numpy NPZ)
   :class: legacy-block

   The legacy pipeline ships two NPZ siblings together: a coefficient file and its source master-dark cube.

   * **Example (coefficients)**: ``masterdarks_coeffs_2025070x_blue.npz``
   * **Example (cube)**: ``masterdarks_2025070x_blue.npz``
   * **Naming**: ``masterdarks_coeffs_YYYYMMxx_[blue|red].npz`` and ``masterdarks_YYYYMMxx_[blue|red].npz``.
   * **Coefficient NPZ arrays**:

     * ``z0``: ``(4400, 4400)`` float64 (intercept per pixel)
     * ``z1``: ``(4400, 4400)`` float64 (slope per pixel)
     * ``logexptime``: ``(7,)`` float64
     * ``ndfilter``: ``(7,)`` float64

   * **Cube NPZ arrays**:

     * ``cube``: ``(4400, 4400, 7)`` float32
     * ``exposuretimes``: ``(7,)`` float64
     * ``ndfilter``: ``(7,)`` float64

   * **Code reference**: ``reduce/recipes/make_coeffs_from_masterdarks.py``

.. topic:: DRAGONS Format (FITS)
   :class: dragons-block

   * **Example**: ``20250707T164838Z_DDDDE_b_0120_darkCoefficients.fits``
   * **Naming**: ``YYYYMMDDTHHmmSSZ_DDDDE_[b|r]_nnnn_darkCoefficients.fits``
   * **Structure**: The ``.data`` extension is a ``(1, 1)`` placeholder;
     the coefficients live in ``.COEFF_Z0``, ``.COEFF_Z1``, and the ``.LOGEXPTIME`` table.

      .. code-block:: pycon

         >>> ad = astrodata.open('20250707T164838Z_DDDDE_b_0120_darkCoefficients.fits')
         >>> ad.info()
         Filename: 20250707T164838Z_DDDDE_b_0120_darkCoefficients.fits
         Tags: 120s BLUE CAL DARK DARK_COEFF GEMINI MAROONX NORTH OVERSCAN_SUBTRACTED
             PREPARED PROCESSED SPECT

         Pixels Extensions
         Index  Content                  Type              Dimensions     Format
         [ 0]   science                  NDAstroData       (1, 1)         float32
                   .variance             ADVarianceUncerta (4400, 4400)   float32
                   .mask                 ndarray           (4400, 4400)   uint16
                   .COEFF_Z0             ndarray           (4400, 4400)   float64
                   .COEFF_Z1             ndarray           (4400, 4400)   float64
                   .LOGEXPTIME           Table             (7, 3)         n/a

         Other Extensions
                        Type        Dimensions
         .EXPOSUREMETER Table       (35, 3)
         .HISTORY       Table       (6, 4)
         .PROVENANCE    Table       (4, 4)

   * **Processing recipe**: ``recipes_DARK.py::makeDarkCoefficients``

Synthetic Darks
---------------

.. topic:: Legacy Format (FITS, synthetic)
   :class: legacy-block

   * **Example**: ``2025070x_masterdark_mean_DDDDE_b_0300.fits`` (300 s blue)
   * **Naming**: ``YYYYMMxx_masterdark_mean_DDDDE_[b|r]_nnnn.fits``.
   * **Structure**: single-extension ``PrimaryHDU`` with a ``(4400, 4400)`` float32 array.

.. topic:: DRAGONS Format (FITS)
   :class: dragons-block

   * **Example**: ``20250717T144308Z_SOOOE_b_0300_synth_dark.fits``
   * **Naming**: takes the timestamp of the source SOOOE frame and appends ``_synth_dark``.
   * **Structure**: single-arm file. ``(4400, 4400)`` float32.

      .. code-block:: pycon

         >>> ad = astrodata.open('20250717T144308Z_SOOOE_b_0300_synth_dark.fits')
         >>> ad.info()
         Filename: 20250717T144308Z_SOOOE_b_0300_synth_dark.fits
         Tags: 300s BLUE CAL DARK DARK_SYNTH GEMINI MAROONX NORTH PREPARED PROCESSED
             SPECT

         Pixels Extensions
         Index  Content                  Type              Dimensions     Format
         [ 0]   science                  NDAstroData       (4400, 4400)   float32
                   .variance             ADVarianceUncerta (4400, 4400)   float32

         Other Extensions
                        Type        Dimensions
         .EXPOSUREMETER Table       (1027, 3)
         .HISTORY       Table       (2, 4)
         .PROVENANCE    Table       (2, 4)

   * **Tag distinguisher**: ``DARK_SYNTH``.
   * **Processing recipe**: ``recipes_ECHELLE_SPECT.py::makeSyntheticDark``
   * **Storage**: written under the ``processed_dark`` caltype in the calibration store.

Flat Field Products
===================

Raw Flats
---------

.. topic:: Legacy Format (FITS)
   :class: legacy-block

   * **Example**: ``20250701T172509Z_DFFFD_b_0007.fits``
   * **Naming**: ``YYYYMMDDTHHmmSSZ_SSSSS_[b|r]_nnnn.fits``

     * ``SSSSS``: fiber configuration (``DFFFD`` for science fibers 2, 3, 4; ``DDDDF`` for fiber 5 only)
     * ``nnnn``: exposure time in integer seconds

   * **Structure**: single-extension ``PrimaryHDU`` with a ``(4400, 4400)`` ``uint16`` array.

.. topic:: DRAGONS Format (FITS/MEF bundle)
   :class: dragons-block

   * **Example**: ``N20250701M6126.fits``
   * **Naming**: ``NYYYYMMDDMiiii.fits``
   * **Structure**: bundle format containing both arms.

      .. code-block:: pycon

         >>> import astrodata
         >>> import maroonx_instruments
         >>> ad = astrodata.open('N20250701M6126.fits')
         >>> ad.info()
         Filename: N20250701M6126.fits
         Tags: BUNDLE CAL FLAT GEMINI MAROONX NORTH RAW SPECT UNPREPARED

         Pixels Extensions
         Index  Content                  Type              Dimensions     Format
         [ 0]   science                  NDAstroData       (4400, 4400)   uint16
         [ 1]   science                  NDAstroData       (4400, 4400)   uint16

         Other Extensions
                        Type        Dimensions
         .EXPOSUREMETER Table       (250, 3)

   * **Processing**: use the ``processBundle`` recipe to separate the two arms before reduction.


Master Flats
------------

.. topic:: Legacy Format (HDF5)
   :class: legacy-block

   * **Example**: ``20250701T17_masterflat_backgroundsubtracted_FFFFF_b_0007.hdf``
   * **Naming**: ``YYYYMMDDTHH_masterflat_backgroundsubtracted_FFFFF_[b|r]_nnnn.hdf``
   * **Root attributes**: ``flat``, ``rawfile`` (paths to source files).
   * **File structure**:

     .. code-block:: text

        stripe_indices/
          fiber                                    # Pixel-to-fiber map
          order                                    # Pixel-to-order map
        box_extraction/fiber_N/order_M/            # 1D box-summed spectra
        optimal_extraction/fiber_N/order_M/        # Optimally extracted spectra
        optimal_var/fiber_N/order_M/               # Variance arrays
        wavelengths/fiber_N/order_M/               # Wavelength grids
        wavelengths_static/fiber_N/order_M/        # Static wavelength grids
        extraction_parameters/fiber_N/order_M/     # Extraction params
        header/                                    # 284 FITS header keys stored as HDF5 attrs

   * **Code references**: ``reduce/recipes/make_master_flats.py``, ``reduce/recipes/backgroundfit.py``, ``reduce/extraction.py``

.. topic:: DRAGONS Format (FITS)
   :class: dragons-block

   * **Example**: ``20250701T172509Z_DDDDF_b_0007_DFFFF_flat.fits``
   * **Naming**: takes the filename of the first raw flat in the input list and appends ``_DFFFF_flat``
   * **Structure**: single-arm file, ``OVERSCAN_SUBTRACTED`` and ``OVERSCAN_TRIMMED`` to ``(4072, 3954)``. 

      .. code-block:: pycon

         >>> ad = astrodata.open('20250701T172509Z_DDDDF_b_0007_DFFFF_flat.fits')
         >>> ad.info()
         Filename: 20250701T172509Z_DDDDF_b_0007_DFFFF_flat.fits
         Tags: 7s BLUE CAL FLAT GEMINI MAROONX NORTH OVERSCAN_SUBTRACTED OVERSCAN_TRIMMED
             PREPARED PROCESSED SPECT

         Pixels Extensions
         Index  Content                  Type              Dimensions     Format
         [ 0]   science                  NDAstroData       (4072, 3954)   float32
                   .variance             ADVarianceUncerta (4072, 3954)   float32
                   .mask                 ndarray           (4072, 3954)   uint16
                   .INDEX_FIBER          ndarray           (4072, 3954)   int64
                   .INDEX_ORDER          ndarray           (4072, 3954)   int64
                   .REMOVED_STRIPES      ndarray           (4, 6)         float64
                   .STRIPES_FIBERS       ndarray           (4,)           int64
                   .STRIPES_ID           Table             (24, 34)       n/a

         Other Extensions
                        Type        Dimensions
         .EXPOSUREMETER Table       (13, 3)
         .HISTORY       Table       (3, 4)
         .PROVENANCE    Table       (2, 4)

   * **Processing recipe**: ``recipes_FLAT_SPECT.py::makeProcessedFlat``


Wavelength Calibration Products
===============================

Raw Wavelength Calibrations
---------------------------

.. topic:: Legacy Format (FITS)
   :class: legacy-block

   * **Example**: ``20250717T163124Z_DEEEE_b_0010.fits``
   * **Naming**: ``YYYYMMDDTHHmmSSZ_SSSSS_[b|r]_nnnn.fits``, one file per arm.

     * ``SSSSS`` = ``DEEEE`` (etalon), ``DTTTT``/``DTTTE``/``DTTTD`` (ThAr),
       or ``DLLLL``/``DLLLE``/``DLLLD``/``DEEEL`` (LFC)

   * **Structure**: single-extension ``PrimaryHDU`` with a ``(4400, 4400)`` ``uint16`` array.

.. topic:: DRAGONS Format (FITS/MEF bundle)
   :class: dragons-block

   * **Example**: ``N20250717M5948.fits`` (etalon, ``DEEEE``)
   * **Naming**: ``NYYYYMMDDMiiii.fits`` (same bundle format as raw darks, flats, and science).
   * **Structure**: bundle format containing both arms.

      .. code-block:: pycon

         >>> import astrodata
         >>> import maroonx_instruments
         >>> ad = astrodata.open('N20250717M5948.fits')
         >>> ad.info()
         Filename: N20250717M5948.fits
         Tags: BUNDLE CAL ETALON GEMINI MAROONX NORTH RAW SPECT UNPREPARED WAVECAL

         Pixels Extensions
         Index  Content                  Type              Dimensions     Format
         [ 0]   science                  NDAstroData       (4400, 4400)   uint16
         [ 1]   science                  NDAstroData       (4400, 4400)   uint16

         Other Extensions
                        Type        Dimensions
         .EXPOSUREMETER Table       (257, 3)


Static Wavelength Solutions
---------------------------

.. topic:: Legacy Format (HDF5)
   :class: legacy-block

   The static solution is stored inside each per-date DEEEE HDF5, evaluated from the reference
   peak model at ``wl_combined_final_etalon_peakmodel_2020.hdf``.

   * **Example**: ``20250709T020701Z_DEEEE_b_0010.hdf``
   * **File structure**:

     .. code-block:: text

        wavelengths_static/fiber_N/order_M/          # N in {2, 3, 4, 5}, M in {91..124} (blue)
          # each order is a (3954,) float64 array of wavelengths [nm]

.. topic:: DRAGONS Format (FITS)
   :class: dragons-block

   In DRAGONS the reference peak model lives as a lookup file
   (``lookups/WLS/REFWAVELENGTH_[b|r].fits``), and ``staticWavelengthSolution()`` evaluates it
   into ``WLS_STATIC_FIBER_*`` extensions on the wavecal file.

   * **Example**: ``20250717T163124Z_DEEEE_b_0010_wavecal.fits``
   * **Structure**: ``(n_orders, n_samples)`` per fiber. Blue arm is ``(34, 3954)``, red arm is
     ``(28, 4036)``. Fiber 1 stays empty as ``(1, 1)``.

      .. code-block:: pycon

         >>> ad.info()
         ...
         [ 0]   science                  NDAstroData       (4072, 3954)   float32
                   ...
                   .WLS_STATIC_FIBER_1   ndarray           (1, 1)         float64
                   .WLS_STATIC_FIBER_2   ndarray           (34, 3954)     float64
                   .WLS_STATIC_FIBER_3   ndarray           (34, 3954)     float64
                   .WLS_STATIC_FIBER_4   ndarray           (34, 3954)     float64
                   .WLS_STATIC_FIBER_5   ndarray           (34, 3954)     float64


Dynamic Wavelength Solutions
----------------------------

.. topic:: Legacy Format (HDF5)
   :class: legacy-block

   * **Example**: ``20250709T020701Z_DEEEE_b_0010.hdf``
   * **File structure**:

     .. code-block:: text

        box_extraction/fiber_N/                      # per-order box-summed spectra
        etalon_peak_parameters/
          peaks                                      # Gaussian fit params
          polynomials                                # polynomial coefficients
        extracted_stripes/fiber_N/                   # rectified 2D stripes
        extraction_parameters/fiber_N/
        peak_data/                                   # pandas DataFrame, multi-index (fiber, order, m)
        stripe_indices/
          fiber                                      # pixel->fiber map
          order                                      # pixel->order map
        wavelengths/fiber_N/order_M/                 # soft link -> wavelengths_dynamic/
        wavelengths_dynamic/fiber_N/order_M/
        wavelengths_static/fiber_N/order_M/
        header/                                      # 253 FITS keys as HDF5 attrs

   * **Instrument drift keywords** (on ``/header`` attrs, stored as byte-strings with units):

     .. code-block:: text

        Instrument_Drift = b'632.4 m/s'
        Drift_Fiber2     = b'623.6 m/s'
        Drift_Fiber3     = b'627.9 m/s'
        Drift_Fiber4     = b'629.8 m/s'

   * **Code reference**: ``analyze/recipes/batch_etalon_spline_wls.py``

.. topic:: DRAGONS Format (FITS)
   :class: dragons-block

   * **Example**: ``20250717T163124Z_DEEEE_b_0010_wavecal.fits``
   * **Naming**: ``YYYYMMDDTHHmmSSZ_DEEEE_[b|r]_nnnn_wavecal.fits``
   * **Structure**: single-arm file, ``OVERSCAN_SUBTRACTED`` and ``OVERSCAN_TRIMMED`` to
     ``(4072, 3954)``. Wavelength solutions, extracted spectra, peak fits, and stripe indices
     all live on extension 0.

      .. code-block:: pycon

         >>> ad = astrodata.open('20250717T163124Z_DEEEE_b_0010_wavecal.fits')
         >>> ad.info()
         Filename: 20250717T163124Z_DEEEE_b_0010_wavecal.fits
         Tags: 10s BLUE CAL ETALON GEMINI MAROONX NORTH OVERSCAN_SUBTRACTED
             OVERSCAN_TRIMMED PREPARED PROCESSED SPECT WAVECAL

         Pixels Extensions
         Index  Content                  Type              Dimensions     Format
         [ 0]   science                  NDAstroData       (4072, 3954)   float32
                   .variance             ADVarianceUncerta (4072, 3954)   float32
                   .mask                 ndarray           (4072, 3954)   uint16
                   .BOX_REDUCED_FIBER_1  ndarray           (1, 1)         float64
                   .BOX_REDUCED_FIBER_2  ndarray           (34, 3954)     float64
                   .BOX_REDUCED_FIBER_3  ndarray           (34, 3954)     float64
                   .BOX_REDUCED_FIBER_4  ndarray           (34, 3954)     float64
                   .BOX_REDUCED_FIBER_5  ndarray           (34, 3954)     float64
                   .BOX_REDUCED_FLAT_1   ndarray           (1, 1)         float64
                   .BOX_REDUCED_FLAT_2   ndarray           (34, 3954)     float64
                   .BOX_REDUCED_FLAT_3   ndarray           (34, 3954)     float64
                   .BOX_REDUCED_FLAT_4   ndarray           (34, 3954)     float64
                   .BOX_REDUCED_FLAT_5   ndarray           (34, 3954)     float64
                   .BOX_REDUCED_VAR_1    ndarray           (1, 1)         float64
                   .BOX_REDUCED_VAR_2    ndarray           (34, 3954)     float64
                   .BOX_REDUCED_VAR_3    ndarray           (34, 3954)     float64
                   .BOX_REDUCED_VAR_4    ndarray           (34, 3954)     float64
                   .BOX_REDUCED_VAR_5    ndarray           (34, 3954)     float64
                   .BPM_FIBER_1          ndarray           (1, 1)         float64
                   .BPM_FIBER_2          ndarray           (34, 3954)     int64
                   .BPM_FIBER_3          ndarray           (34, 3954)     int64
                   .BPM_FIBER_4          ndarray           (34, 3954)     int64
                   .BPM_FIBER_5          ndarray           (34, 3954)     int64
                   .INDEX_FIBER          ndarray           (4072, 3954)   int64
                   .INDEX_ORDER          ndarray           (4072, 3954)   int64
                   .PEAKS                Table             (69753, 11)    n/a
                   .PEAK_DATA            Table             (69653, 16)    n/a
                   .POLY                 Table             (136, 10)      n/a
                   .REDUCED_ORDERS_FIBER ndarray           (1, 1)         float64
                   .REDUCED_ORDERS_FIBER ndarray           (34,)          float64
                   .REDUCED_ORDERS_FIBER ndarray           (34,)          float64
                   .REDUCED_ORDERS_FIBER ndarray           (34,)          float64
                   .REDUCED_ORDERS_FIBER ndarray           (34,)          float64
                   .STRIPES_ID           Table             (24, 34)       n/a
                   .WLS_DYNAMIC_FIBER_1  ndarray           (1, 1)         float64
                   .WLS_DYNAMIC_FIBER_2  ndarray           (34, 3954)     float64
                   .WLS_DYNAMIC_FIBER_3  ndarray           (34, 3954)     float64
                   .WLS_DYNAMIC_FIBER_4  ndarray           (34, 3954)     float64
                   .WLS_DYNAMIC_FIBER_5  ndarray           (34, 3954)     float64
                   .WLS_STATIC_FIBER_1   ndarray           (1, 1)         float64
                   .WLS_STATIC_FIBER_2   ndarray           (34, 3954)     float64
                   .WLS_STATIC_FIBER_3   ndarray           (34, 3954)     float64
                   .WLS_STATIC_FIBER_4   ndarray           (34, 3954)     float64
                   .WLS_STATIC_FIBER_5   ndarray           (34, 3954)     float64

         Other Extensions
                        Type        Dimensions
         .EXPOSUREMETER Table       (257, 3)
         .HISTORY       Table       (3, 4)
         .PROVENANCE    Table       (3, 4)

   * **Processing recipe**: ``recipes_DYNAMIC_WAVECAL.py::makeDynamicWavecal``


Simultaneous Wavelength Solutions
---------------------------------

.. topic:: Legacy Format (HDF5)
   :class: legacy-block

   Applied inside a science-frame HDF5 (``SOOOE`` fiber setup), fiber 5 is the reference etalon
   and fiber 6 stores the etalon wavelengths after drift correction.

   * **Example**: ``20250817T131207Z_SOOOE_b_0300.hdf``
   * **Added groups** (on top of the dynamic-wavecal tree shown above):

     .. code-block:: text

        wavelengths_simultaneous/fiber_N/order_M/    # N in {2, 3, 4, 5, 6}
        wavelengths/                                 # soft link -> wavelengths_simultaneous/

   * **Instrument drift keywords** on the science HDF5:

     .. code-block:: text

        Instrument_Drift = b'789.8 m/s'
        Relative_Drift   = b'-7.0 m/s'

   * **Code reference**: ``analyze/recipes/batch_science_spline_wls_dynamic.py``

.. topic:: DRAGONS Format (FITS)
   :class: dragons-block

   Applied to science frames by the ``reduce`` recipe, producing ``WLS_SIMULTANEOUS_FIBER_*``
   extensions. Fiber 1 stays empty; fibers 2, 3, 4 carry the drift-corrected science
   wavelengths; fiber 5 holds the reference etalon.

   * **Example**: ``20241124T041907Z_SOOOE_r_0300_reduced.fits``
   * **Naming**: ``YYYYMMDDTHHmmSSZ_SOOOE_[b|r]_nnnn_reduced.fits``

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

     Extension names are truncated by column width in ``ad.info()``. In full they are
     ``WLS_SIMULTANEOUS_FIBER_1`` through ``WLS_SIMULTANEOUS_FIBER_5``.

   * **Processing recipe**: ``recipes_ECHELLE_SPECT.py::reduce``

Science Frame Products
======================

Raw Science Frames
------------------

.. topic:: Legacy Format (FITS)
   :class: legacy-block

   * **Example**: ``20250717T144308Z_SOOOE_b_0300.fits``
   * **Naming**: ``YYYYMMDDTHHmmSSZ_SOOOE_[b|r]_nnnn.fits``, one file per arm.

     * ``SOOOE``: fiber configuration (Sky, Object, Object, Object, Etalon)
     * ``[b|r]``: blue or red camera arm
     * ``nnnn``: exposure time in integer seconds

   * **Structure**: single-extension ``PrimaryHDU`` with a ``(4400, 4400)`` ``uint16`` array.

.. topic:: DRAGONS Format (FITS/MEF bundle)
   :class: dragons-block

   * **Example**: ``N20250717M5299.fits``
   * **Naming**: ``NYYYYMMDDMiiii.fits`` (same bundle format as raw darks, flats, and wavecals).
   * **Structure**: MEF bundle with both arms.

      .. code-block:: pycon

         >>> import astrodata
         >>> import maroonx_instruments
         >>> ad = astrodata.open('N20250717M5299.fits')
         >>> ad.info()
         Filename: N20250717M5299.fits
         Tags: 300s BUNDLE GEMINI MAROONX NORTH RAW SCI SPECT UNPREPARED

         Pixels Extensions
         Index  Content                  Type              Dimensions     Format
         [ 0]   science                  NDAstroData       (4400, 4400)   uint16
         [ 1]   science                  NDAstroData       (4400, 4400)   uint16

         Other Extensions
                        Type        Dimensions
         .EXPOSUREMETER Table       (1027, 3)

         >>> ad.fiber_setup()
         ['Sky', 'Target', 'Target', 'Target', 'Etalon']

   * **Processing**: use the ``processBundle`` recipe to separate the two arms before reduction.


Reduced Science Frames
----------------------

.. topic:: Legacy Format (Pandas HDF5 store)
   :class: legacy-block

   The final reduced science product is a Pandas ``HDFStore`` (opened with ``pd.HDFStore``), not a
   plain HDF5 tree. Both arms live in one file.

   * **Example**: ``20250717T144308Z_SOOOE_x_0300.hd5``
   * **Naming**: ``YYYYMMDDTHHmmSSZ_SOOOE_x_nnnn.hd5``. The ``x`` marks a combined-arm file, and the
     extension is ``.hd5`` (not ``.hdf`` as used elsewhere in the legacy pipeline).
   * **File structure**:

     .. code-block:: text

        /spec_[blue|red]/                # DataFrame, MultiIndex (Fiber, Order)
                                         # 170 rows = 5 fibers x 34 orders (blue arm)
          columns: box_extraction, wavelengths, optimal_extraction, optimal_var
          each cell: (3954,) float64 array (blue arm)

        /blaze_[blue|red]/               # DataFrame (empty placeholder)
        /header_[blue|red]/              # Series, 313 FITS keys per arm

   * **Code reference**: ``analyze/recipes/combine_science_fibers.py`` (adds the combined science
     fiber to the per-arm HDF5 intermediates before this final store is written).

.. topic:: DRAGONS Format (FITS)
   :class: dragons-block

   * **Example**: ``20250717T144308Z_SOOOE_b_0300_reduced.fits``
   * **Naming**: ``YYYYMMDDTHHmmSSZ_SOOOE_[b|r]_nnnn_reduced.fits``, one file per arm.
   * **Structure**: single-arm file, ``OVERSCAN_SUBTRACTED`` and ``OVERSCAN_TRIMMED`` to
     ``(4072, 3954)``. Box and optimal extractions, wavelength solutions, per-fiber BPMs, order
     lists, and etalon fit tables all live on extension 0.

      .. code-block:: pycon

         >>> ad = astrodata.open('20250717T144308Z_SOOOE_b_0300_reduced.fits')
         >>> ad.info()
         Filename: 20250717T144308Z_SOOOE_b_0300_reduced.fits
         Tags: 300s BLUE GEMINI MAROONX NORTH OVERSCAN_SUBTRACTED OVERSCAN_TRIMMED
             PREPARED PROCESSED PROCESSED_SCIENCE SCI SPECT

         Pixels Extensions
         Index  Content                  Type              Dimensions     Format
         [ 0]   science                  NDAstroData       (4072, 3954)   float32
                   .variance             ADVarianceUncerta (4072, 3954)   float32
                   .mask                 ndarray           (4072, 3954)   uint16
                   .BOX_REDUCED_FIBER_1  ndarray           (1, 1)         float64
                   .BOX_REDUCED_FIBER_2  ndarray           (34, 3954)     float64
                   .BOX_REDUCED_FIBER_3  ndarray           (34, 3954)     float64
                   .BOX_REDUCED_FIBER_4  ndarray           (34, 3954)     float64
                   .BOX_REDUCED_FIBER_5  ndarray           (34, 3954)     float64
                   .BOX_REDUCED_FLAT_1   ndarray           (1, 1)         float64
                   .BOX_REDUCED_FLAT_2   ndarray           (34, 3954)     float64
                   .BOX_REDUCED_FLAT_3   ndarray           (34, 3954)     float64
                   .BOX_REDUCED_FLAT_4   ndarray           (34, 3954)     float64
                   .BOX_REDUCED_FLAT_5   ndarray           (34, 3954)     float64
                   .BOX_REDUCED_VAR_1    ndarray           (1, 1)         float64
                   .BOX_REDUCED_VAR_2    ndarray           (34, 3954)     float64
                   .BOX_REDUCED_VAR_3    ndarray           (34, 3954)     float64
                   .BOX_REDUCED_VAR_4    ndarray           (34, 3954)     float64
                   .BOX_REDUCED_VAR_5    ndarray           (34, 3954)     float64
                   .BPM_FIBER_1          ndarray           (1, 1)         float64
                   .BPM_FIBER_2          ndarray           (34, 3954)     int64
                   .BPM_FIBER_3          ndarray           (34, 3954)     int64
                   .BPM_FIBER_4          ndarray           (34, 3954)     int64
                   .BPM_FIBER_5          ndarray           (34, 3954)     int64
                   .INDEX_FIBER          ndarray           (4072, 3954)   int64
                   .INDEX_ORDER          ndarray           (4072, 3954)   int64
                   .OPTIMAL_REDUCED_FIBE ndarray           (1, 1)         float64
                   .OPTIMAL_REDUCED_FIBE ndarray           (34, 3954)     float64
                   .OPTIMAL_REDUCED_FIBE ndarray           (34, 3954)     float64
                   .OPTIMAL_REDUCED_FIBE ndarray           (34, 3954)     float64
                   .OPTIMAL_REDUCED_FIBE ndarray           (1, 1)         float64
                   .OPTIMAL_REDUCED_FIBE ndarray           (34, 3954)     float64
                   .OPTIMAL_REDUCED_VAR_ ndarray           (1, 1)         float64
                   .OPTIMAL_REDUCED_VAR_ ndarray           (34, 3954)     float64
                   .OPTIMAL_REDUCED_VAR_ ndarray           (34, 3954)     float64
                   .OPTIMAL_REDUCED_VAR_ ndarray           (34, 3954)     float64
                   .OPTIMAL_REDUCED_VAR_ ndarray           (1, 1)         float64
                   .OPTIMAL_REDUCED_VAR_ ndarray           (34, 3954)     float64
                   .PEAKS                Table             (17438, 11)    n/a
                   .POLY                 Table             (34, 10)       n/a
                   .REDUCED_ORDERS_FIBER ndarray           (1, 1)         float64
                   .REDUCED_ORDERS_FIBER ndarray           (34,)          float64
                   .REDUCED_ORDERS_FIBER ndarray           (34,)          float64
                   .REDUCED_ORDERS_FIBER ndarray           (34,)          float64
                   .REDUCED_ORDERS_FIBER ndarray           (34,)          float64
                   .REDUCED_ORDERS_FIBER ndarray           (34,)          float64
                   .STRIPES_ID           Table             (24, 34)       n/a
                   .WLS_SIMULTANEOUS_FIB ndarray           (1, 1)         float64
                   .WLS_SIMULTANEOUS_FIB ndarray           (34, 3954)     float64
                   .WLS_SIMULTANEOUS_FIB ndarray           (34, 3954)     float64
                   .WLS_SIMULTANEOUS_FIB ndarray           (34, 3954)     float64
                   .WLS_SIMULTANEOUS_FIB ndarray           (34, 3954)     float64
                   .WLS_SIMULTANEOUS_FIB ndarray           (34, 3954)     float64
                   .WLS_STATIC_FIBER_1   ndarray           (1, 1)         float64
                   .WLS_STATIC_FIBER_2   ndarray           (34, 3954)     float64
                   .WLS_STATIC_FIBER_3   ndarray           (34, 3954)     float64
                   .WLS_STATIC_FIBER_4   ndarray           (34, 3954)     float64
                   .WLS_STATIC_FIBER_5   ndarray           (34, 3954)     float64

         Other Extensions
                        Type        Dimensions
         .EXPOSUREMETER Table       (1027, 3)
         .HISTORY       Table       (4, 4)
         .PROVENANCE    Table       (6, 4)

     Extension names are truncated by column width in ``ad.info()``. In full,
     ``OPTIMAL_REDUCED_FIBE...`` covers ``OPTIMAL_REDUCED_FIBER_1`` through
     ``OPTIMAL_REDUCED_FIBER_6``, and the same six-fiber pattern applies to
     ``OPTIMAL_REDUCED_VAR_*`` and ``WLS_SIMULTANEOUS_FIBER_*``.

   * **Fiber roles**: fibers 2, 3, 4 are the science targets; fiber 5 is the simultaneous etalon
     reference; fiber 6 is the combined science spectrum (2+3+4); fiber 1 (sky) stays as a
     ``(1, 1)`` placeholder.
   * **Processing recipe**: ``recipes_ECHELLE_SPECT.py::reduce``


Extraction Products Comparison
------------------------------

Compact mapping from legacy Pandas HDF5 keys to DRAGONS FITS extensions:

.. list-table:: Science Frame Data Product Mapping
   :header-rows: 1
   :widths: 45 45

   * - Legacy Pandas HDF5 (``.hd5``)
     - DRAGONS FITS extension
   * - ``/spec_[blue|red]`` col ``box_extraction``
     - ``.BOX_REDUCED_FIBER_N``
   * - ``/spec_[blue|red]`` col ``optimal_extraction``
     - ``.OPTIMAL_REDUCED_FIBER_N``
   * - ``/spec_[blue|red]`` col ``optimal_var``
     - ``.OPTIMAL_REDUCED_VAR_N``
   * - ``/spec_[blue|red]`` col ``wavelengths``
     - ``.WLS_SIMULTANEOUS_FIBER_N``
   * - ``/header_[blue|red]`` Series
     - PHU header keys



Configuration Files
===================

DRAGONS moves the standing calibration content from the legacy configuration
files into typed lookup files under ``maroonxdr/maroonx/lookups/``. Two legacy
sources feed this:

* ``config_[b|r].hdf`` (one pair per arm, updated occasionally): bad pixel map,
  stripe identification, and the evaluated static wavelength grid.
* ``wl_combined_final_etalon_peakmodel_2020.hdf`` (shared reference peak model):
  the polynomial fit that the evaluated static grid is derived from.

Each subsection below maps one legacy chunk to its DRAGONS lookup.

Bad Pixel Masks
---------------

.. topic:: Legacy Format (HDF5)
   :class: legacy-block

   Bad pixel masks live inside the shared per-arm configuration files. There is a
   single ``config_[b|r].hdf`` pair per arm at any given time, updated
   occasionally when the calibration changes.

   * **Location**: ``.../MaroonX_spectra_reduced/Maroonx_configfiles/YYYYMMxx/config_[b|r].hdf``
   * **BPM-related datasets**:

     .. code-block:: text

        bad_pixel_map           # (4400, 4400) int64 (big-endian)
        valid                   # (4400, 4400) int64  valid-region mask

   The same ``config_[b|r].hdf`` file also carries the stripe-identification
   content and the reference wavelength grids, covered in the next two
   subsections.

.. topic:: DRAGONS Format (FITS)
   :class: dragons-block

   Bad pixel masks are standalone FITS files under
   ``maroonxdr/maroonx/lookups/BPM/``:

   .. code-block:: text

      BPM_[b|r]_0000.fits   # PrimaryHDU (4400, 4400) int64

   .. code-block:: pycon

      >>> from astropy.io import fits
      >>> hdul = fits.open('maroonxdr/maroonx/lookups/BPM/BPM_b_0000.fits')
      >>> hdul.info()
      Filename: BPM_b_0000.fits
      No.    Name      Ver    Type      Cards   Dimensions   Format
        0  PRIMARY       1 PrimaryHDU       8   (4400, 4400)   int64

   * **Registered in**: ``maroonxdr/maroonx/lookups/maskdb.py``


Stripe Identification
---------------------

.. topic:: Legacy Format (HDF5)
   :class: legacy-block

   Stripe metadata is stored inside the same ``config_[b|r].hdf`` file. The
   groups ``find_stripes/``, ``identify_stripes/``, and
   ``correct_image_orientation/`` are empty as datasets but carry their
   parameters as HDF5 **attributes**:

   .. code-block:: text

      find_stripes/                          # attrs: parameters for stripe finding
        deg_polynomial       = 5
        gauss_filter_sigma   = 3.5
        min_peak             = 0.008

      identify_stripes/                      # attrs: known stripe positions
        positions            = (170, 3) int array   # rows = (fiber, order, y-position on detector)
                                                    # 170 = 5 fibers x 34 orders
                                                    # first row: [1, 100, 1341]

      correct_image_orientation/             # attrs: detector orientation flags
        fliplr               = True
        flipud               = True

      used_thar_lines/                       # 33 datasets, one per echelle order
                                             # each: (N_lines, 2) float64  (pixel, wavelength_nm)

.. topic:: DRAGONS Format (FITS)
   :class: dragons-block

   Files under ``maroonxdr/maroonx/lookups/SID/``:

   .. code-block:: text

      SID_[b|r].fits    
      README

   .. code-block:: pycon

      >>> from astropy.io import fits
      >>> hdul = fits.open('maroonxdr/maroonx/lookups/SID/SID_b.fits')
      >>> hdul.info()
      Filename: SID_b.fits
      No.    Name      Ver    Type      Cards   Dimensions   Format
        0  PRIMARY       1 PrimaryHDU       9   ()
        1  SID           1 BinTableHDU     15   170R x 3C   [I, I, I]

   Columns of the ``SID`` extension:

   * ``identify_fiber`` (int16): fiber number (1-5)
   * ``fiber_order`` (int16): echelle order number
   * ``fiber_position`` (int16): y-position on detector, pixels

   The ``SID`` table is a direct FITS-native port of the legacy
   ``identify_stripes/positions`` attribute. 
   The ``find_stripes/`` params live as ``findStripes``
   defaults in ``parameters_maroonx_2D.py``.


Reference Wavelengths
---------------------

See also **Static Wavelength Solutions** in the Wavelength Calibration
Products section above.

.. topic:: Legacy Format (HDF5)
   :class: legacy-block

   Two related HDF5 sources.

   * **Reference peak model**: ``wl_combined_final_etalon_peakmodel_2020.hdf``.

     .. code-block:: text

        dispersion/
          parameter                          # JSON-encoded lmfit Parameters

        wls_blue/fiber_N/                    # N in {2, 3, 4}
          maxx                               # scalar int64 = 3954
          orders                             # int64
          wavelengths                        # float64 [nm]
          weights                            # float64
          x_norm                             # float64
          poly_deg_x                         # scalar int64 = 16
          poly_deg_y                         # scalar int64 = 23

        wls_red/fiber_N/                     # same schema, maxx=4036

   * **Evaluated static grids** inside the same ``config_[b|r].hdf`` file:

     .. code-block:: text

        wavelengths_static/fiber_N/order_M/  # N in {1, 2, 3, 4, 5}, M in {91..124} (blue)
                                             # each order is a (3954,) float64 wavelength grid [nm]


.. topic:: DRAGONS Format (FITS)
   :class: dragons-block

   Two lookup files under ``maroonxdr/maroonx/lookups/WLS/``.

   * **Reference peak model**: ``REFWAVELENGTH_[b|r].fits`` mirrors the legacy
     ``_peakmodel_2020.hdf``.

     .. code-block:: text

        PARAMETERS   Table (99, 8)           # polynomial fit params
                                             # columns: Name, Value, Min, Max, Stderr, Vary, Expr, Brute_Step
        FIBER_N      Table (1, 4)            # N in {2, 3, 4}
                                             # columns: XNORM, ORDERS, WEIGHTS, WAVELEN

   * **Pre-evaluated static solution**: ``WLSTAT_[b|r].fits``.

     .. code-block:: pycon

        >>> from astropy.io import fits
        >>> hdul = fits.open('maroonxdr/maroonx/lookups/WLS/WLSTAT_b.fits')
        >>> hdul.info()
        Filename: WLSTAT_b.fits
        No.    Name      Ver    Type      Cards   Dimensions   Format
          0  PRIMARY       1 PrimaryHDU       9   ()
          1  FIBER_1       1 BinTableHDU     77   3954R x 34C   [D, D, D, D, ...]
          2  FIBER_2       1 BinTableHDU     77   3954R x 34C   [D, D, D, D, ...]
          3  FIBER_3       1 BinTableHDU     77   3954R x 34C   [D, D, D, D, ...]
          4  FIBER_4       1 BinTableHDU     77   3954R x 34C   [D, D, D, D, ...]
          5  FIBER_5       1 BinTableHDU     77   3954R x 34C   [D, D, D, D, ...]

     Each ``FIBER_N`` has 34 columns named by echelle order number (``100``,
     ``101``, ..., ``133``), each holding 3954 float64 wavelength values in nm.

   ``WLSTAT_[b|r].fits`` is the ready-to-use per-fiber static wavelength
   solution and is what ``staticWavelengthSolution()`` reads at runtime.
   ``REFWAVELENGTH_[b|r].fits`` holds the underlying peak model in case
   re-evaluation is ever needed.

   * **Code reference**: ``utils/ref_wls_hdf2fits.py``


