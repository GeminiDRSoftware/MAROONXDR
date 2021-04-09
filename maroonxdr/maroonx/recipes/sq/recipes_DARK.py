"""
Recipes available to data with tags ['MAROONX', 'CAL', 'DARK'].
Default is "makeProcessedDark"
"""

recipe_tags = set(['MAROONX', 'CAL', 'DARK'])

def makeProcessedDark(p):
    """
    This recipe performs the standardization and corrections needed to convert
    the raw input dark images into a single master dark image. This output
    processed dark is stored on disk using storeProcessedDark and has a name
    equal to the name of the first input bias image with "_dark.fits" appended.

    Parameters
    ----------
    p : PrimitivesCORE object
        A primitive set matching the recipe_tags.
    """
    #p.somestuff()

    p.checkArm()
    p.checkND()
    p.prepare()
    # p.addDQ(add_illum_mask=False)  # need to get MX BPM
    # p.overscanCorrect()  # I think this almost works but we have horizontal sections for the overscan
    # (e.g. rawdata[:2200,2050:2200] is an overscan region being applied to rawdata[:2040,:2040])

    p.addVAR(read_noise=True)
    p.ADUToElectrons()
    p.addVAR(poisson_noise=True)
    p.stackDarks()
    # want to do (each pixel) 2 sigma clipped mean between frames
    # (i.e. mean_frame, _, _ = astropy.stats.sigma_clipped_stats(data_cube,axis=2,sigma=2.0)
    p.storeProcessedDark()

    return

_default = makeProcessedDark


