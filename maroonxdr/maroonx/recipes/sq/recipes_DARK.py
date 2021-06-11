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
    The background in an un-illuminated frame is very low for exposure times of
    less than 900s and likely doesn't warrant a dark subtraction, however
    most science frames are taken with the simultaneous calibration fiber
    illuminated with the FP etalon. The extended wings of the etalon reach
    into one of the science fibers with a few 10 counts. To remove these wings
    and any broad diffuse (illumination independent) background, DDDDE frames
    are taken in daytime to construct a DDDDE master dark. These darks are specific
    for different exposure times (or ND filter settings) and should be taken
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
    # p.correctImageOrientation()
    # p.addDQ(add_illum_mask=False)  # need to get MX BPM
    p.addVAR(read_noise=True)
    p.overscanCorrect()  # I think this almost works but we have horizontal sections for the overscan
    # (e.g. rawdata[:2200,2050:2200] is an overscan region being applied to rawdata[:2040,:2040])
    p.ADUToElectrons()
    p.addVAR(poisson_noise=True)
    p.stackDarks()
    # want to do (each pixel) 2 sigma clipped mean between frames
    # (i.e. mean_frame, _, _ = astropy.stats.sigma_clipped_stats(data_cube,axis=2,sigma=2.0)
    # and save the min value for each pixel across the original frames
    # (i.e. min_frame = np.min(data_cube,axis=2)
    # then identify and replace incorrectly clipped data through the following algorithm
    # diff = mean_frame - min_frame
    # bad_locs = np.where((diff > 2*np.nanmedian(mean_frame)) & (diff_frame > np.abs(0.3*min)) & (mean_frame > 10) & (diff_frame > 10) )
    # if np.shape(bad_locs)[1] > 0:
    #     mean[bad_locs] = min[bad_locs]

    # algorithm is: The difference between sigma-clipped mean and minimum must be larger than 2x the global median
    # and larger than 30% of the minimum. Both the mean and the difference of those pixels must also be greater than
    # 10 counts. This corrects for coincidentally overlapping cosmic ray hits i the input files
    p.storeProcessedDark()

    return

_default = makeProcessedDark


