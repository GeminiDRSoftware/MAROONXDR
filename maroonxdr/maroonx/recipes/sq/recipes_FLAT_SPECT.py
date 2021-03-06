"""
Recipes available to data with tags ['MAROONX', 'CAL', 'FLAT'].
Default is "makeProcessedFlat"
"""

recipe_tags = set(['MAROONX', 'CAL', 'FLAT'])

def makeProcessedFlat(p):
    """
    This recipe performs the standardization and corrections needed to convert
    the raw input flat images into a single stacked flat image. This output
    processed flat is stored on disk using storeProcessedFlat and has a name
    equal to the name of the first input bias image with "_flat.fits" appended.
    A master flatfield is required to perform optimal flux extraction and to
    determine the blaze function. Due to the stability of the spectrograph, one master
    flatfield is typically 'valid' for at least a two-week period. Cross-comparison
    between different master flatfields taken months apart have not been conducted so far.
    Parameters
    ----------
    p : PrimitivesCORE object
        A primitive set matching the recipe_tags.
    """
    p.prepare()
    p.checkArm()
    p.checkND()
    # p.correctImageOrientation()
    # p.addDQ()
    # p.addVAR(read_noise=False)
    # p.overscanCorrect()  # I think this almost works but we have horizontal sections for the overscan
    #     # (e.g. rawdata[:2200,2050:2200] is an overscan region being applied to rawdata[:2040,:2040])

    p.separateFlatStreams()  # creates 'DFFFD_flats' stream and leaves FDDDF flats in main stream

    p.stackFlats()
    p.stackFlats(stream='DFFFD_flats')

    # correct_image_orientation, find_stripes, identify_stripes,
    # need to implement illuminated fiber order tracing here (to mask them for background fitting)

    # need to do diffuse background subtraction here, based on 2D fit of masked frame

    # need to combine DFFFD and FDDDF frames (i.e. just make with masterflat image = np.max([DFFFD_b,FDDDF_b],axis=0)
    p.combineFlatStreams(stream='main', source='DFFFD_flats')
    p.clearStream(stream='DFFFD_flats')
    # run the 5-illuminated-fiber frame through extraction to create a reduced masterflat
    p.storeProcessedFlat()
    return

_default = makeProcessedFlat


