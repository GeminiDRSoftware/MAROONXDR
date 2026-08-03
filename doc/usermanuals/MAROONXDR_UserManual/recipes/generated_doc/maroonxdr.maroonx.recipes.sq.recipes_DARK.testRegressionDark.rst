testRegressionDark
==================

| **Recipe Library**: maroonxdr.maroonx.recipes.sq.recipes_DARK
| **Astrodata Tags**: {'MAROONX', 'CAL', 'DARK'}

Produce a processed dark for regression comparison against legacy.

Test-support recipe. Follows the same processing as makeProcessedDark,
without overscan trimming nor image orientation correction to match the
legacy dark production, and stores the result with a "_regressionDark"
suffix for use in regression tests.

::

    Parameters
    ----------
    p : Primitives object
        A primitive set matching the recipe_tags.

::

    def testRegressionDark(p):
        p.prepare()
        p.checkArm()
        p.checkND()
        p.addDQ()
        p.subtractOverscan()
        # No trim overscan and no correct image orientation
        p.addVAR(read_noise=True, poisson_noise=True)
        p.stackDarks(scale_mode='first_frame', lsigma=2.0, hsigma=2.0)
        p.storeProcessedDark(suffix='_regressionDark')

