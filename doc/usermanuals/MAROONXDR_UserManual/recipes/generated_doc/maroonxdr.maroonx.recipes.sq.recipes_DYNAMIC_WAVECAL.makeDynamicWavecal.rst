makeDynamicWavecal
==================

| **Recipe Library**: maroonxdr.maroonx.recipes.sq.recipes_DYNAMIC_WAVECAL
| **Astrodata Tags**: {'MAROONX', 'WAVECAL'}

Process MAROON-X 2D etalon exposures into a dynamic wavelength solution.

This is done in the following steps:

1. Utilizing the relevant flat, the fibers and orders are traced and
   identified from the frame. These are extracted into sparse arrays.
2. Using box extraction, these are converted into 1D spectra. Box
   extraction is the simple summation of all spatial pixels in a given
   fiber and order combination. The trace of the 'box' is taken from
   the processed flat field.
3. The extracted 1D etalon lines are fitted to determine their
   centroids. This involves identifying the peaks and fitting them to a
   box convolved with two Gaussians. The locations of the peak centers
   are stored. The width of the peaks and the Gaussian sigmas vary very
   slowly along an order and are modeled by low-order polynomials.
4. The dynamic wavelength solution is computed by fitting the 1D pixel
   positions and the wavelengths of the etalon lines with a 30 knot
   cubic spline. Pixel positions are identified by comparison to the
   static wavelength solution, loaded from a lookup file, which is
   accurate to about 500 meters per second. This measures the drift of
   the spectrograph with time and restores an accuracy of 10 to 20
   centimeters per second.

The result is stored in the calibration database as a processed wavecal.

::

    Parameters
    ----------
    p : Primitives object
        A primitive set matching the recipe_tags.

::

    def makeDynamicWavecal(p):
        p.prepare()
        p.checkArm()
        p.addDQ()  # just placeholder until MX is in caldb

        p.subtractOverscan()
        p.trimOverscan()
        p.correctImageOrientation()
        p.addVAR(read_noise=True, poisson_noise=True)

        p.extractStripes()
        p.boxExtraction()
        p.getPeaksAndPolynomials()

        p.staticWavelengthSolution()
        p.fitAndApplyEtalonWls()
        p.storeProcessedWavecal(suffix='_wavecal')

