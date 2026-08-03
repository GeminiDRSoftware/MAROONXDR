makeBlaze
=========

| **Recipe Library**: maroonxdr.maroonx.recipes.sq.recipes_FLAT_SPECT
| **Astrodata Tags**: {'MAROONX', 'CAL', 'FLAT'}

Measure the blaze function for each fiber of a processed masterflat.

This recipe fits a spline to the box-extracted flat spectrum of each
fiber to model the blaze curve, then normalises each order so the peak
equals 1. The result is stored as a BLAZE_FIBER_N extension for each
fiber N and written to disk with a "_blaze" suffix.

::

    Parameters
    ----------
    p : Primitives object
        A primitive set matching the recipe_tags.

::

    def makeBlaze(p):
        p.checkMaster()
        p.extractStripes()
        # p.boxExtraction()
        p.optimalExtraction()
        p.measureBlaze(n_knots=50)
        p.writeOutputs(suffix='_blaze')

