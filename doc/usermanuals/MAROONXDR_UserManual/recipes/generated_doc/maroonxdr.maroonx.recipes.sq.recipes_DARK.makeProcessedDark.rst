makeProcessedDark
=================

| **Recipe Library**: maroonxdr.maroonx.recipes.sq.recipes_DARK
| **Astrodata Tags**: {'DARK', 'MAROONX', 'CAL'}

Convert raw MAROON-X dark frames into a single processed dark.

The input dark frames are checked for arm and ND filter consistency,
then stacked into one processed dark image. The result is stored on
disk by storeProcessedDark under the name of the first input dark
with "_dark.fits" appended.

The background in an un-illuminated frame is very low for exposure
times below 900s and would hardly warrant a dark subtraction on its
own. However, most science frames are taken with the simultaneous
calibration fiber illuminated by the FP etalon, whose extended line
wings reach into one of the science fibers at the level of a few tens
of counts. To remove these wings, together with any broad diffuse
(illumination independent) background, DDDDE frames are taken in the
daytime and stacked into a DDDDE processed dark.

::

    Parameters
    ----------
    p : Primitives object
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

