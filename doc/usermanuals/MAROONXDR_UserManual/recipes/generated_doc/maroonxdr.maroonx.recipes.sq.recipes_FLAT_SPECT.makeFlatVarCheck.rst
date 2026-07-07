makeFlatVarCheck
================

| **Recipe Library**: maroonxdr.maroonx.recipes.sq.recipes_FLAT_SPECT
| **Astrodata Tags**: {'CAL', 'FLAT', 'MAROONX'}

Check if variance extensions are correctly computed on stacked flats.

This recipe does not find, identify, or define any stripes. It also does
not remove stray light. Mostly used to test if variance is being computed
correctly for a stack of images.

::

    Parameters
    ----------
    p : PrimitivesCORE object
        A primitive set matching the recipe_tags.

    Returns
    -------
    Creates test frames with variance added.

::

    def makeFlatVarCheck(p):
        p.prepare()
        p.checkArm()
        # p.checkND()
        p.addDQ()
        p.subtractOverscan()
        p.trimOverscan()
        p.correctImageOrientation()
        p.addVAR(read_noise=True, poisson_noise=True)
        # Creates 'DFFFD_flats' stream and leaves FDDDF flats in main stream
        p.separateFlatStreams()
        p.stackFlats(stream='DFFFD_flats', suffix='DFFFD')
        p.storeProcessedFlat(stream='DFFFD_flats', suffix='_varAddedStack')

