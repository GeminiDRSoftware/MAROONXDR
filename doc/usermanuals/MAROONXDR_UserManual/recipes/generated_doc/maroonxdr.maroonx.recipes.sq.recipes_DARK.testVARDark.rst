testVARDark
===========

| **Recipe Library**: maroonxdr.maroonx.recipes.sq.recipes_DARK
| **Astrodata Tags**: {'CAL', 'MAROONX', 'DARK'}

Produce a dark frame with an additional variance plane.

Can be used to add a variance plane to a singular dark frame. The default
recipe adds variance planes to all the dark frames while stacking, and also
outputs a dark frame with a variance plane added to it.

::

    Parameters
    ----------
    p : PrimitivesCORE object
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

