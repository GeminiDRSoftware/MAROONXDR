"""
Recipes available to data with tags ['MAROONX', 'CAL', 'FLAT'].

Default is "makeProcessedFlat".
TODO: Add a step to get the 1D Spectra for the fibers to use in Wavecal and
Science.
"""

recipe_tags = {'MAROONX', 'CAL', 'FLAT'}

def makeProcessedFlat(p):
    """
    Perform standardization and corrections to convert raw flats to processed.

    This recipe converts the raw input flat images into a single stacked flat
    image. This output processed flat is stored on disk using storeProcessedFlat
    and has a name equal to the name of the first input bias image with
    "_FFFFF_flat.fits" appended. A processed flatfield is required to perform
    optimal flux extraction and to determine the blaze function. Due to the
    stability of the spectrograph, one processed flatfield is typically 'valid'
    for at least a two-week period. Cross-comparison between different processed
    flatfields taken months apart have not been conducted so far.

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
    p.addVAR(read_noise=True,poisson_noise=True)
    # Creates 'DFFFD_flats' stream and leaves FDDDF flats in main stream
    p.separateFlatStreams()
    p.stackFlats(suffix='FDDDF_flats')
    p.stackFlats(stream='DFFFD_flats',suffix='DFFFD')
    # Define stripe info to ultimately remove stray light in each stream
    p.findStripes()
    p.findStripes(stream='DFFFD_flats')
    # Identify stripes based on MX architecture files
    p.identifyStripes(selected_fibers=[1, 5])
    p.identifyStripes(stream='DFFFD_flats', selected_fibers=[2, 3, 4])
    # Defines pixel inclusion for each flat region based on stripe ids
    p.defineFlatStripes()
    p.defineFlatStripes(stream='DFFFD_flats')
    # Remove straylight (requires 2 partial illumination flat sets)
    p.removeStrayLight(stream='main', filter_size=19, box_size=20)
    p.removeStrayLight(stream='DFFFD_flats', filter_size=19, box_size=20)

    # Legacy patch for removeStrayLight
    # p.removeStrayLight_legacyPatch(stream='main', filter_size=19, box_size=20)
    # p.removeStrayLight_legacyPatch(stream='DFFFD_flats', filter_size=19, box_size=20) 

    # Combine straylight-removed images
    p.combineFlatStreams(stream='main', stream_2='DFFFD_flats')

    # Remove second stream
    p.clearStream(stream='DFFFD_flats')
    # Re-run find/identify/define routine on combined frame
    p.findStripes()
    p.identifyStripes(selected_fibers=[1, 2, 3, 4, 5])
    p.defineFlatStripes(extract=True)
    
    # Perform optimal extraction on flat field to create 1D spectra
    p.extractStripes()
    p.optimalExtraction(optimal_extraction_fibers=[2, 3, 4, 5])

    p.storeProcessedFlat(suffix='_FFFFF_flat')

_default = makeProcessedFlat

def makeStrayLightCheck(p):
    """
    Check the stray light subtraction in normal flat frame processing.

    Run the straylight_test_prep.py file to generate these.

    Parameters
    ----------
    p : PrimitivesCORE object
        A primitive set matching the recipe_tags.

    Returns
    -------
    Creates test frames with straylight difference and flux at levels just
    before straylight removal. Unit test will independently perform straylight
    removal and compare results.
    """
    p.prepare()
    p.checkArm()
    p.checkND()
    p.addDQ()
    p.subtractOverscan()
    p.trimOverscan()
    p.correctImageOrientation()
    p.addVAR(read_noise=True,poisson_noise=True)
    # Creates 'DFFFD_flats' stream and leaves FDDDF flats in main stream
    p.separateFlatStreams()
    p.stackFlats(suffix='FDDDF_flats')
    p.stackFlats(stream='DFFFD_flats', suffix='DFFFD')
    # Define stripe info to ultimately remove stray light in each stream
    p.findStripes()
    p.findStripes(stream='DFFFD_flats')
    # Identify stripes based on MX architecture files
    p.identifyStripes(selected_fibers=[1, 5])
    p.identifyStripes(stream='DFFFD_flats', selected_fibers=[2, 3, 4])
    # Defines pixel inclusion for each flat region based on stripe ids
    p.defineFlatStripes()
    p.defineFlatStripes(stream='DFFFD_flats')
    # Remove straylight (requires 2 partial illumination flat sets)
    p.removeStrayLight(snapshot=True, filter_size=19, box_size=20)
    p.removeStrayLight(
        stream='DFFFD_flats', snapshot=True, filter_size=19, box_size=20
    )
    p.storeProcessedFlat(stream='DFFFD_flats', suffix='_straylight_flat')
    p.storeProcessedFlat(suffix='_straylight_flat')

def makeFlatVarCheck(p):
    """
    Check if variance extensions are correctly computed on stacked flats.

    This recipe does not find, identify, or define any stripes. It also does
    not remove stray light. Mostly used to test if variance is being computed
    correctly for a stack of images.

    Parameters
    ----------
    p : PrimitivesCORE object
        A primitive set matching the recipe_tags.

    Returns
    -------
    Creates test frames with variance added.
    """
    p.prepare()
    p.checkArm()
    p.checkND()
    p.addDQ()
    p.subtractOverscan()
    p.trimOverscan()
    p.correctImageOrientation()
    p.addVAR(read_noise=True,poisson_noise=True)
    # Creates 'DFFFD_flats' stream and leaves FDDDF flats in main stream
    p.separateFlatStreams()
    p.stackFlats(stream='DFFFD_flats', suffix='DFFFD')
    p.storeProcessedFlat(stream='DFFFD_flats', suffix='_varAddedStack')


def makeProcessedFlatDFFFF(p):
    """
    Perform standardization and corrections to convert raw flats to processed.

    This recipe converts the raw input flat images into a single stacked flat
    image. This output processed flat is stored on disk using storeProcessedFlat
    and has a name equal to the name of the first input bias image with
    "_FFFFF_flat.fits" appended. A processed flatfield is required to perform
    optimal flux extraction and to determine the blaze function. Due to the
    stability of the spectrograph, one processed flatfield is typically 'valid'
    for at least a two-week period. Cross-comparison between different processed
    flatfields taken months apart have not been conducted so far.

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
    # p.trimOverscan()  # noqa: ERA001
    # p.correctImageOrientation()  # noqa: ERA001
    p.addVAR(read_noise=True, poisson_noise=True)
    # Creates 'DFFFD_flats' stream and leaves FDDDF flats in main stream
    p.separateFlatStreams()

    p.stackFlats(stream='main', scale_mode='mean_frame', suffix='DDDDF')
    p.stackFlats(stream='DFFFD_flats', scale_mode='mean_frame', suffix='DFFFD')

    # Subtract overscan is run again in backgroundfit.py on the stacked flats
    p.subtractOverscan(stream='main')
    p.subtractOverscan(stream='DFFFD_flats')

    p.trimOverscan(stream='main')
    p.trimOverscan(stream='DFFFD_flats')

    p.correctImageOrientation(stream='main')
    p.correctImageOrientation(stream='DFFFD_flats')
    # ================================================

    # Define stripe info to ultimately remove stray light in each stream
    p.findStripes(stream='main')
    p.findStripes(stream='DFFFD_flats')

    # Identify stripes based on MX architecture files
    p.identifyStripes(stream='main', selected_fibers=[5])
    p.identifyStripes(stream='DFFFD_flats', selected_fibers=[2, 3, 4])

    # Defines pixel inclusion for each flat region based on stripe ids
    p.defineFlatStripes(stream='main')
    p.defineFlatStripes(stream='DFFFD_flats')

    # Remove straylight (requires 2 partial illumination flat sets)
    p.removeStrayLight(stream='main', filter_size=19, box_size=20)
    p.removeStrayLight(stream='DFFFD_flats', filter_size=19, box_size=20)

    # Legacy patch for removeStrayLight
    # p.removeStrayLight_legacyPatch(stream='main', filter_size=19, box_size=20)
    # p.removeStrayLight_legacyPatch(stream='DFFFD_flats', filter_size=19, box_size=20)

    # Combine straylight-removed images
    p.combineFlatStreams(stream='main', stream_2='DFFFD_flats')
    # Remove second stream
    p.clearStream(stream='DFFFD_flats')

    # Re-run find/identify/define routine on combined frame
    p.findStripes()
    p.identifyStripes(selected_fibers=[2,3,4,5])
    p.defineFlatStripes(extract=True)
    
    # Perform optimal extraction on flat field to create 1D spectra
    # TODO: optimal extraction requires a flat calibration which is the product
    # of this recipe. This needs to be a new recipe that is run after this one (?).
    # p.extractStripes()
    # p.optimalExtraction(optimal_extraction_fibers=[2, 3, 4, 5])

    p.storeProcessedFlat(suffix='_DFFFF_flat')