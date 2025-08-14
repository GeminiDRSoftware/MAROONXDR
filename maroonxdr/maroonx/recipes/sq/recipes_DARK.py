"""
Recipes available to data with tags ['MAROONX', 'CAL', 'DARK'].
Default is "makeProcessedDark"
"""

recipe_tags = set(['MAROONX', 'CAL', 'DARK'])

def makeProcessedDark(p):
    """
    This recipe performs the standardization and corrections needed to convert
    the raw input dark images into a single processed dark image. This output
    processed dark is stored on disk using storeProcessedDark and has a name
    equal to the name of the first input bias image with "_dark.fits" appended.
    The background in an un-illuminated frame is very low for exposure times of
    less than 900s and likely doesn't warrant a dark subtraction, however
    most science frames are taken with the simultaneous calibration fiber
    illuminated with the FP etalon. The extended wings of the etalon reach
    into one of the science fibers with a few 10 counts. To remove these wings
    and any broad diffuse (illumination independent) background, DDDDE frames
    are taken in daytime to construct a DDDDE processed dark. These darks are specific
    for different exposure times (i.e. ND filter settings) and should be taken
    close in time (within a day or two) to the science frame as the etalon
    source brightness can be time variable.

    Parameters
    ----------
    p : PrimitivesCORE object
        A primitive set matching the recipe_tags.
    """
    p.prepare()
    p.checkArm()
    p.checkND()
    p.addDQ()
    p.subtractOverscan()
    # legacy dark coefficients dont trim nor flip
    p.trimOverscan()
    p.correctImageOrientation()

    p.addVAR(read_noise=True, poisson_noise=True)
    p.stackDarks(scale_mode='first_frame', lsigma=2.0, hsigma=2.0)
    p.storeProcessedDark()
    return

_default = makeProcessedDark

def makeDarkCoefficients(p):
    """
    This recipe produces the coefficient arrays 'z0' and 'z1' from a 
    pixel-by-pixel fit of the relationship F = z1 + z0 * log10(Texp). The input 
    is a list of MASTER darks created from raw dark files. The input should 
    cover a wide range of exposure times and can have multiple master dark 
    files with the same exposure time. Master darks with non-standard ND filter 
    settings should be avoided (they will not fit the same relationship as 
    the other frames).
    """
    p.prepare()
    p.checkArm()
    p.checkMaster()
    # p.checkND()
    p.fitDarkCoefficients()
    p.storeProcessedDark(suffix="_darkCoefficients")


def testVARDark(p):
    """
    This recipe produces a dark frame with an additional variance plane added to it.
    Can be used to add a variance plane to a singular dark frame.  The default recipe
    adds variance planes to all the dark frames while stacking, and also outputs a dark 
    frame with a variance plane added to it. 

    Parameters
    ----------
    p : PrimitivesCORE object
        A primitive set matching the recipe_tags.
    """
    p.prepare()
    p.checkArm()
    p.checkND()
    p.addDQ()
    p.subtractOverscan()
    p.trimOverscan()
    p.correctImageOrientation()
    p.addVAR(read_noise=True, poisson_noise=True)
    p.storeProcessedDark(suffix='_varAdded')


def testRegressionDark(p):
    """
    This recipe is used to test the dark frame of the legacy pipeline can be reproduced.
    It is a simplified version of the makeProcessedDark recipe, it does not
    include overscan trimming nor image orientation correction.

    Parameters
    ----------
    p : PrimitivesCORE object
        A primitive set matching the recipe_tags.
    """

    p.prepare()
    p.checkArm()
    p.checkND()
    p.addDQ()
    p.subtractOverscan()
    # No trim overscan and no correct image orientation
    p.addVAR(read_noise=True, poisson_noise=True)
    p.stackDarks(scale_mode='first_frame', lsigma=2.0, hsigma=2.0)
    p.storeProcessedDark(suffix='_regressionDark')
    return