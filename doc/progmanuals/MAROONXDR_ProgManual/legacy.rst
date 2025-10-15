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

**GOA Archive Format** (January 2025 onwards):

  ``NYYYYMMDDMnnnn.fits``. Example: ``N20241115M0001.fits``

Where ``N`` = Gemini North, ``M`` = MAROON-X, and ``nnnn`` is a sequential counter
within each UT day. This is a MEF container for the red and blue files, which can be restored
to the native format using the ``processBundle`` recipe.


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
   * **Tags**: ``{'GEMINI', 'BUNDLE', 'DARK', 'UNPREPARED', 'RAW', 'MAROONX'}``

Master Darks
------------

.. topic:: Legacy Format (FITS)
   :class: legacy-block

   * **Naming**: ``YYYYMMDDTHH_masterdark_mean_DDDDE_[b|r]_nnnn.fits``

   * **Structure**: Single ImageHDU with 4400×4400 pixels (float64)

   * **Standard Exposure Times**: 60, 120, 300, 600, 900, 1200, 1800 seconds
   * **Code Reference**: ``reduce/recipes/make_master_darks.py``

.. topic:: DRAGONS Format (FITS/MEF)
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
         Tags:

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

   * **Tags**: ``{'GEMINI', '[BLUE|RED]', 'DARK', 'CAL', 'PROCESSED', 'MAROONX'}``

Dark Coefficient Files
----------------------

.. topic:: Legacy Format (NumPy NPZ)
   :class: legacy-block

   Coefficient files enable synthetic dark generation as a function of exposure time:

   * **Naming**: ``masterdarks_coeffs_YYYYMMxx_[blue|red].npz``
   * **Format**: NumPy compressed archive (NPZ) with three arrays:

     * ``z0``: Slope coefficients (4400×4400, float64)
     * ``z1``: Intercept coefficients (4400×4400, float64)
     * ``logexptime``: log(exposure times) used for fitting

   * **Code Reference**: ``reduce/recipes/make_coeffs_from_masterdarks.py``

.. topic:: DRAGONS Format (FITS)
   :class: dragons-block

   Coefficient files store the same log-linear model using FITS table extensions:

   * **Naming**: ``YYYYMMDDTHHmmSSZ_DDDDE_[b|r]_nnnn_darkCoefficients.fits``
   * **Structure**: Single-arm with coefficient extensions:

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


Flat Field Products
===================

.. todo::
   Document raw flats, master flats, and straylight-corrected flats

Wavelength Calibration Products
================================

.. todo::
   Document wavelength solution files and comparison with legacy format

Extracted Spectra Products
===========================

.. todo::
   Document 1D extracted spectra format and comparison with legacy

Order Tracing Products
======================

.. todo::
   Document order trace solutions and comparison with legacy

Bad Pixel Masks
===============

.. todo::
   Document BPM format and comparison with legacy

Science Frame Products
======================

.. todo::
   Document reduced science frames and comparison with legacy

File Naming Conventions
========================

.. todo::
   Compare legacy and DRAGONS file naming schemes

Metadata and Provenance
========================

.. todo::
   Discuss differences in metadata storage and provenance tracking