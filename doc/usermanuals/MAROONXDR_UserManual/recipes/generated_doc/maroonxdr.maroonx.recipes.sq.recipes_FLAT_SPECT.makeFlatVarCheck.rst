makeFlatVarCheck
================

| **Recipe Library**: maroonxdr.maroonx.recipes.sq.recipes_FLAT_SPECT
| **Astrodata Tags**: {'CAL', 'MAROONX', 'FLAT'}

Check if variance extensions are correctly computed on stacked flats.

Test-support recipe. The DFFFD flats are stacked and stored with a
"_varAddedStack" suffix, without any stripe definition or stray light
removal. Mostly used to test if variance is being computed correctly
for a stack of images.

::

    Parameters
    ----------
    p : Primitives object
        A primitive set matching the recipe_tags.

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

