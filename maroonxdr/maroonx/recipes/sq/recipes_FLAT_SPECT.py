"""
Recipes available to data with tags ['MAROONX', 'CAL', 'FLAT'].
Default is "makeProcessedFlat"
TODO: Add a step to get the 1D Spectra for the fibers to use in Wavecal and Science
"""

recipe_tags = set(['MAROONX', 'CAL', 'FLAT'])

def makeProcessedFlat(p):
    """
    This recipe performs the standardization and corrections needed to convert
    the raw input flat images into a single stacked flat image. This output
    processed flat is stored on disk using storeProcessedFlat and has a name
    equal to the name of the first input bias image with "_FFFFF_flat.fits" appended.
    A processed flatfield is required to perform optimal flux extraction and to
    determine the blaze function. Due to the stability of the spectrograph, one processed
    flatfield is typically 'valid' for at least a two-week period. Cross-comparison
    between different processed flatfields taken months apart have not been conducted so far.
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
    p.separateFlatStreams()  # creates 'DFFFD_flats' stream and leaves FDDDF flats in main stream
    p.stackFlats(suffix='FDDDF_flats')
    p.stackFlats(stream='DFFFD_flats',suffix='DFFFD')
    p.findStripes()  # define stripe info to ultimately remove stray light in each stream
    p.findStripes(stream='DFFFD_flats')
    p.identifyStripes(selected_fibers='1,0,0,0,5') # identify stripes based on MX architecture files
    p.identifyStripes(stream='DFFFD_flats',selected_fibers='0,2,3,4,0')
    p.defineFlatStripes()  # defines pixel inclusion for each flat region based on stripe ids
    p.defineFlatStripes(stream='DFFFD_flats')
    p.removeStrayLight(filter_size = 21, box_size = 21)  # remove straylight from frame (this is why 2 partial illumination flat sets are necessary)
    p.removeStrayLight(stream='DFFFD_flats', filter_size = 21, box_size = 21)
    p.combineFlatStreams(stream='main', source='DFFFD_flats')  # combine straylight-removed images
    p.clearStream(stream='DFFFD_flats') # remove second stream
    p.findStripes()  # re-run find/identify/define routine on combined frame
    p.identifyStripes(selected_fibers='1,2,3,4,5')
    p.defineFlatStripes(extract=True)
    # run the 5-illuminated-fiber frame through extraction to create a reduced processed flat
    p.storeProcessedFlat(suffix='_FFFFF_flat')
    return

_default = makeProcessedFlat

def makeStrayLightCheck(p):
    """
    This recipe is utilized to check the stray light subtraction that is made
    in the normal processing of a series of flat frames.  Run the straylight_test_prep.py
    file to generate these.
    Parameters
    ----------
    p : PrimitivesCORE object
        A primitive set matching the recipe_tags.
    Returns
    -------
    creates test frames with straylight difference and flux at levels just before straylight removal
    unit test will independently preform straylight removal and compare results
    """
    p.prepare()
    p.checkArm()
    p.checkND()
    p.addDQ()
    p.overscanCorrect()
    p.correctImageOrientation()
    p.addVAR(read_noise=True,poisson_noise=True)
    p.separateFlatStreams()  # creates 'DFFFD_flats' stream and leaves FDDDF flats in main stream
    p.stackFlats(suffix='FDDDF_flats')
    p.stackFlats(stream='DFFFD_flats', suffix='DFFFD')
    p.findStripes()  # define stripe info to ultimately remove stray light in each stream
    p.findStripes(stream='DFFFD_flats')
    p.identifyStripes(selected_fibers='1,0,0,0,5')  # identify stripes based on MX architecture files
    p.identifyStripes(stream='DFFFD_flats', selected_fibers='0,2,3,4,0')
    p.defineFlatStripes()  # defines pixel inclusion for each flat region based on stripe ids
    p.defineFlatStripes(stream='DFFFD_flats')
    p.removeStrayLight(snapshot=True, filter_size = 21, box_size = 21)  # remove straylight from frame (this is why 2 partial illumination flat sets are necessary)
    p.removeStrayLight(stream='DFFFD_flats', snapshot=True, filter_size = 21, box_size = 21)
    p.storeProcessedFlat(stream='DFFFD_flats', suffix='_straylight_flat')
    p.storeProcessedFlat(suffix='_straylight_flat')

def makeFlatVarCheck(p):
    """
    This recipe is used to check if the variance extensions are correctly being 
    computed to a stack of flat images.  It does not find, idntify, or define any stripes.
    It also does not remove stray light.  Mostly used to test if variance is being 
    computed correctly for a stack of images.
    Parameters
    ----------
    p : PrimitivesCORE object
        A primitive set matching the recipe_tags.
    Returns
    -------
    creates test frames with variance added
    """
    p.prepare()
    p.checkArm()
    p.checkND()
    p.addDQ()
    p.overscanCorrect()
    p.correctImageOrientation()
    p.addVAR(read_noise=True,poisson_noise=True)
    p.separateFlatStreams()  # creates 'DFFFD_flats' stream and leaves FDDDF flats in main stream
    p.stackFlats(stream='DFFFD_flats', suffix='DFFFD')
    p.storeProcessedFlat(stream='DFFFD_flats', suffix='_varAddedStack')