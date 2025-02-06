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
    p.overscanCorrect()
    p.correctImageOrientation()
    p.addVAR(read_noise=True,poisson_noise=True)
    p.stackDarks()  # see parameters file for sig-clip choices and overwrites
    p.storeProcessedDark()

    return

_default = makeProcessedDark

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
    p.overscanCorrect()
    p.correctImageOrientation()
    p.addVAR(read_noise=True,poisson_noise=True)
    p.storeProcessedDark(suffix='_varAdded')
