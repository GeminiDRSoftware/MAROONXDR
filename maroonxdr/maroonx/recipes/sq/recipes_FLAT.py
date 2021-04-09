"""
Recipes available to data with tags ['MAROONX', 'CAL', 'BIAS'].
Default is "makeProcessedBias"
"""

recipe_tags = set(['MAROONX', 'CAL', 'BIAS'])

def makeProcessedFLAT(p):
    """
    This recipe performs the standardization and corrections needed to convert
    the raw input flat images into a single stacked flat image. This output
    processed flat is stored on disk using storeProcessedFlat and has a name
    equal to the name of the first input bias image with "_flat.fits" appended.

    Parameters
    ----------
    p : PrimitivesCORE object
        A primitive set matching the recipe_tags.
    """

    # need to do analysis of DFFFD and FDDDF frames in parallel (as possible for either arm)
    p.checkArm()
    p.checkND()
    p.prepare()
    # p.addDQ()
    p.addVAR(read_noise=True)
    # p.overscanCorrect()
    p.stackFlats()

    # need to implement illuminated fiber order tracing here (to mask them later)

    # need to do background subtraction here, based on 2D fit of masked frame

    # need to combine DFFFD and FDDDF frames (i.e. just make with np.max([DFFFD_b,FDDDF_b],axis=0))

    p.storeProcessedFlat()
    return

_default = makeProcessedFLAT


