makeProcessedDark
=================

| **Recipe Library**: maroonxdr.maroonx.recipes.sq.recipes_DARK
| **Astrodata Tags**: {'DARK', 'CAL', 'MAROONX'}

Perform standardization and corrections to convert raw darks to processed.

This recipe converts the raw input dark images into a single processed dark
image. This output processed dark is stored on disk using storeProcessedDark
and has a name equal to the name of the first input bias image with
"_dark.fits" appended. The background in an un-illuminated frame is very low
for exposure times of less than 900s and likely doesn't warrant a dark
subtraction, however most science frames are taken with the simultaneous
calibration fiber illuminated with the FP etalon. The extended wings of the
etalon reach into one of the science fibers with a few 10 counts. To remove
these wings and any broad diffuse (illumination independent) background,
DDDDE frames are taken in daytime to construct a DDDDE processed dark. These
darks are specific for different exposure times (i.e. ND filter settings)
and should be taken close in time (within a day or two) to the science frame
as the etalon source brightness can be time variable.

::

    Parameters
    ----------
    p : PrimitivesCORE object
        A primitive set matching the recipe_tags.

::

    def makeProcessedDark(p):
        p.prepare()
        p.checkArm()
        p.checkND()
        p.addDQ()
        p.subtractOverscan()
        # legacy dark coefficients dont trim nor flip
        # p.trimOverscan()  # noqa: ERA001
        # p.correctImageOrientation()  # noqa: ERA001

        p.addVAR(read_noise=True, poisson_noise=True)
        p.stackDarks(scale_mode='first_frame', lsigma=2.0, hsigma=2.0)
        p.storeProcessedDark()

