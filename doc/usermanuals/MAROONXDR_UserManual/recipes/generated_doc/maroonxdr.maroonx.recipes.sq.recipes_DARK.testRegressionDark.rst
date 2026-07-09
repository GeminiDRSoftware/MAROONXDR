testRegressionDark
==================

| **Recipe Library**: maroonxdr.maroonx.recipes.sq.recipes_DARK
| **Astrodata Tags**: {'CAL', 'MAROONX', 'DARK'}

Test that the legacy pipeline dark frame can be reproduced.

This is a simplified version of the makeProcessedDark recipe that does not
include overscan trimming nor image orientation correction.

::

    Parameters
    ----------
    p : PrimitivesCORE object
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

