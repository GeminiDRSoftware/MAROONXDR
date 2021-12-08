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
    # p.addDQ()  # need to get bpm frame read in correctly
    p.overscanCorrect()
    p.correctImageOrientation()
    p.addVAR(read_noise=True,poisson_noise=True)
    p.separateFlatStreams()  # creates 'DFFFD_flats' stream and leaves FDDDF flats in main stream
    p.stackFlats(suffix='FDDDF_flats')
    p.stackFlats(stream='DFFFD_flats',suffix='DFFFD')
    p.findStripes()
    p.findStripes(stream='DFFFD_flats')
    p.identifyStripes(selected_fibers='1,0,0,0,5')
    p.identifyStripes(stream='DFFFD_flats',selected_fibers='0,2,3,4,0')
    p.defineFlatStripes()
    p.defineFlatStripes(stream='DFFFD_flats')
    p.removeStrayLight()
    p.removeStrayLight(stream='DFFFD_flats')
    p.combineFlatStreams(stream='main', source='DFFFD_flats')

    p.clearStream(stream='DFFFD_flats')
    # run the 5-illuminated-fiber frame through extraction to create a reduced masterflat
    p.storeProcessedFlat()
    return

_default = makeProcessedFlat


