.. _maroonxdr_user_instrument:

*******************
Instrument and Data
*******************

General Description
===================

MAROON-X is a fiber-fed, cross-dispersed echelle spectrograph installed at the
Gemini-North Observatory. Its primary science goal is the detection of
Earth-sized planets in the habitable zones of mid to late M dwarf stars
through high-precision radial velocity measurements, with a target precision
of 1 m/s.

The instrument has no movable parts and operates in a single fixed
configuration. Incoming light is split into two arms that are
recorded simultaneously on independent CCD detectors:

- **Blue arm**: 491-670 nm
- **Red arm**: 649-920 nm

The blue detector is a 4400×4400 pixel, 16-bit array read through four
amplifiers (quadrants Q1–Q4), with a gain of 2.72 e⁻/DN and read noise of
1.14 (variance units). The red detector is also 4400×4400 pixels, read through
two amplifiers (R1–R2), with a gain of 2.74 e⁻/DN and read noise of 1.63
(variance units).

Raw data delivered from the Gemini Observatory Archive (GOA) arrives as a
bundle — a single Multi-Extension FITS (MEF) file containing both arms
(``NYYYYMMDDMnnnn.fits``). The first pipeline step is to split this bundle
into separate blue and red arm files for independent processing.

Description of the Modes
=========================

MAROON-X has only one observing mode. At each wavelength, five fibers are
arranged along the cross-dispersion direction on the detector:

- **Fiber 1**: Off-target sky background fiber
- **Fibers 2, 3, 4**: Three pupil-sliced fractions of the on-sky science target
- **Fiber 5**: Simultaneous calibration fiber (fed calibration light during
  science observations)

Each fiber traces a stripe across the detector for each echelle order. The
blue arm covers 34 orders per fiber and the red arm covers 28 orders per fiber.

The fiber illumination pattern for each frame is encoded as a five-character
string — one character per fiber — in the file header and used by the pipeline
for frame classification:

.. list-table::
   :header-rows: 1
   :widths: 10 40

   * - Code
     - Illumination source
   * - ``D``
     - Dark (no illumination)
   * - ``F``
     - Flat field lamp
   * - ``O``
     - Object (stellar target)
   * - ``S``
     - Sky background (fiber 1 only)
   * - ``E``
     - Fabry-Perot etalon comb
   * - ``T``
     - ThAr hollow-cathode arc lamp
   * - ``I``
     - Iodine cell
   * - ``L``
     - Laser Frequency Comb (LFC)

Common frame types and their fiber patterns are summarised in the table below:

.. list-table::
   :header-rows: 1
   :widths: 15 15 40

   * - Frame type
     - Pattern
     - Purpose
   * - Science
     - ``SOOOE``
     - Standard science: sky in fiber 1, target in fibers 2–4, etalon in fiber 5
   * - Dark
     - ``DDDDE``
     - Dark calibration with etalon in simultaneous calibration fiber
   * - Flat (science fibers)
     - ``DFFFD``
     - Flat field for science fibers 2, 3, 4
   * - Flat (outer fibers)
     - ``FDDDF``
     - Flat field for fibers 1 and 5
   * - Etalon wavecal
     - ``DEEEE``
     - Dynamic (nightly) wavelength calibration
   * - ThAr wavecal
     - ``DTTTT``
     - Static wavelength calibration
   * - LFC wavecal
     - ``DLLLL``
     - Alternative wavelength calibration with Laser Frequency Comb

Required Calibration and Associated Observations
================================================

The following calibration frames are required to reduce MAROON-X data:

**Master Darks (DDDDE)**

Dark frames are taken with all science fibers dark and the etalon illuminating
fiber 5. This is necessary because during science observations the etalon wings
contribute a few tens of counts into adjacent fibers. Dark frames must be
matched to the science frame by **exposure time** and by the **ND filter
position** of the simultaneous calibration fiber, and should be taken within
one to two days of the science frames due to variability in the etalon source
brightness.

Standard exposure times are 60, 120, 300, 600, 900, 1200, and 1800 seconds.
Multiple master darks at different exposure times can be combined into
pixel-by-pixel coefficient files (``z0``, ``z1``), parameterising the
relationship ``F = z1 + z0 × log10(Texp)``. These coefficient files allow
synthetic master darks to be generated for any required exposure time.

**Master Flats (DFFFD + FDDDF)**

Two complementary sets of flat frames are required, because no single
illumination pattern can expose all five fibers simultaneously without the
wings of illuminated fibers contaminating adjacent dark fibers for straylight
removal:

- ``DFFFD``: illuminates science fibers 2, 3, 4 only.
- ``FDDDF``: illuminates the outer fibers 1 and 5 only.

Each partial-illumination flat has its 2D straylight background removed before
the two sets are combined into a single five-fiber (``FFFFF``) master flat. The
master flat provides the fiber/order trace polynomials used by all subsequent
extractions and the blaze function for flux normalisation.

Due to the spectrograph's stability, one processed master flat is typically
valid for at least two weeks.

**Static Wavelength Calibration (ThAr, DTTTT)**

ThAr hollow-cathode lamp frames establish an absolute wavelength scale accurate
to approximately 500 m/s. This calibration has been performed only once in the
instrument's lifetime. The resulting solution is stored as precomputed lookup
files (``WLSTAT_b.fits`` and ``WLSTAT_r.fits``) and is not regenerated through
the pipeline.

**Dynamic Wavelength Calibration (Etalon, DEEEE)**

Fabry-Perot etalon frames are taken every night of observing to track
instrumental drift relative to the static wavelength solution. Fitting etalon
peak positions to a cubic spline restores wavelength accuracy from ~500 m/s
(static) to 10-20 cm/s.

**Simultaneous Wavelength Calibration (Etalon in Fiber 5)**

During standard science observations (``SOOOE``), fiber 5 is simultaneously
illuminated by the etalon. This provides a continuous record of instrumental
drift at the time of each science exposure and is used in the
``applyWavelengthSolution`` step to correct the science fibers (2, 3, 4).

Important Instrument Characteristics or Issues
==============================================

**Sky fiber (fiber 1) illumination failure**

The original flat-field design called for two complementary illumination
patterns — ``DFFFD`` (science fibers 2, 3, 4) and ``FDDDF`` (outer fibers 1
and 5) — so that all five fibers could be traced and the combined master flat
would be a fully illuminated ``FFFFF`` frame. In practice, there have been
issues illuminating the sky fiber (fiber 1), and ``FDDDF`` flats have been
replaced by ``DDDDF`` frames that illuminate only fiber 5. The pipeline
therefore accepts ``DDDDF`` as a substitute for ``FDDDF`` when building the
master flat.
