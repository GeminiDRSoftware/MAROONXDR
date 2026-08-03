testVARDark
===========

| **Recipe Library**: maroonxdr.maroonx.recipes.sq.recipes_DARK
| **Astrodata Tags**: {'MAROONX', 'CAL', 'DARK'}

Produce a dark frame with an additional variance plane.

Test-support recipe. Each input dark is processed individually (overscan
subtracted, trimmed, orientation corrected) without stacking, and stored
with a variance plane computed from read noise and Poisson noise. Intended
for checking the variance computation on a single dark frame.

::

    Parameters
    ----------
    p : Primitives object
        A primitive set matching the recipe_tags.

::

    def testVARDark(p):
        p.prepare()
        p.checkArm()
        p.checkND()
        p.addDQ()
        p.subtractOverscan()
        p.trimOverscan()
        p.correctImageOrientation()
        p.addVAR(read_noise=True, poisson_noise=True)
        p.storeProcessedDark(suffix='_varAdded')

