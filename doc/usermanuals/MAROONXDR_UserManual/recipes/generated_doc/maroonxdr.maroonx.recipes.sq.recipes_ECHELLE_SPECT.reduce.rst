reduce
======

| **Recipe Library**: maroonxdr.maroonx.recipes.sq.recipes_ECHELLE_SPECT
| **Astrodata Tags**: {'MAROONX', 'SCI'}

Process MAROON-X science echelle spectrum with tracing and extraction.

This recipe: (1) traces and identifies the fibers and orders using a 2D
processed flat, (2) performs both regular (aka 'box') and optimal
extraction to produce 1D extracted spectra, (3) computes a drift
corrected wavelength solution for the science fibers, and (4) combines
the science fibers and calculates the barycentric velocity correction.
The result is stored on disk with a "_reduced" suffix.

Tracing and identifying fibers and orders is done on a (preferably
background subtracted) 2D processed flat retrieved from the calibration
database. During the stripe extraction a matching processed dark is
subtracted from the science fibers, while the sim cal fiber gets its
straylight removed instead.

Box extraction is the simple summation of all spatial pixels in a given
fiber and order combination. Optimal extraction is by default only
applied to the science fibers 2, 3 and 4.

The wavelength calibration fits the etalon lines of the sim cal fiber
5, loads the static wavelength solution from a lookup file, and applies
a drift corrected solution to the science fibers by comparison with a
processed wavecal etalon frame retrieved from the calibration database.
The barycentric velocity correction is computed from the exposure meter
flux-weighted timestamps and stored in header keywords.

::

    Parameters
    ----------
    p : Primitives object
        A primitive set matching the recipe_tags.

::

    def reduce(p):
        p.prepare()
        p.checkArm()
        p.addDQ()  # just placeholder until MX is in caldb
        p.overscanCorrect()
        p.correctImageOrientation()
        p.addVAR(read_noise=True, poisson_noise=True)
        p.extractStripes(dark_subtraction_skip_fibers=[5], straylight_removal_fibers=[5])
        p.optimalExtraction()
        p.getPeaksAndPolynomials(fibers=(5,))
        p.staticWavelengthSolution()
        p.applyWavelengthSolution(fibers=(2, 3, 4), ref_fiber=5)
        p.combineFibers()
        p.barycentricCorrection()
        p.storeProcessedScience(suffix='_reduced')

