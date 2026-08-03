"""
Recipes available to data with tags ['MAROONX', 'CAL', 'FLAT'].

Default is "makeProcessedFlat".
"""

recipe_tags = {'MAROONX', 'CAL', 'FLAT'}
blocked_tags = {'BUNDLE'}


def makeProcessedFlat(p):
    """
    Convert raw MAROON-X flat frames into a single processed flat.

    The input flats are separated into two streams based on their fiber
    illumination pattern: FDDDF flats (fibers 1 and 5 illuminated) in the
    main stream and DFFFD flats (fibers 2, 3 and 4 illuminated) in a
    second stream. Each stream is stacked, its fiber stripes are traced
    and identified, and its stray light is removed. The two streams are
    then combined into a fully illuminated FFFFF flat by taking the
    pixel-by-pixel maximum, the stripe tracing is re-run on the combined
    frame, and 1D spectra are extracted for the illuminated fibers. The
    result is stored on disk by storeProcessedFlat under the name of the
    first input flat with "_FFFFF_flat.fits" appended.

    A processed flatfield is required to perform optimal flux extraction
    and to determine the blaze function. Due to the stability of the
    spectrograph, one processed flatfield is typically valid for at least
    a two-week period. Cross-comparisons between processed flatfields
    taken months apart have not been conducted so far.

    Parameters
    ----------
    p : Primitives object
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

    Test-support recipe. Mirrors the makeProcessedFlatDFFFF flow up to
    removeStrayLight, run with snapshot option enabled: the SCI data is
    left at the level just before straylight removal and the removed
    straylight is saved as the STRAYLIGHT_DIFFERENCE extension. A unit
    test independently performs the straylight removal and compares its
    result against the sum of the SCI and STRAYLIGHT_DIFFERENCE extensions.

    Parameters
    ----------
    p : Primitives object
        A primitive set matching the recipe_tags.
    """
    p.prepare()
    p.checkArm()
    p.addDQ()
    p.subtractOverscan()
    p.addVAR(read_noise=True, poisson_noise=True)
    # Creates 'DFFFD_flats' stream and leaves DDDDF flats in main stream
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
    p.removeStrayLight(stream='main', snapshot=True, filter_size=19, box_size=20)
    p.removeStrayLight(stream='DFFFD_flats', snapshot=True, filter_size=19, box_size=20)

    p.writeOutputs(stream='DFFFD_flats', suffix='_straylight_flat', strip=True)
    p.writeOutputs(suffix='_straylight_flat', strip=True)


# old recipe - set for deprecation
def makeFlatVarCheck(p):
    """
    Check if variance extensions are correctly computed on stacked flats.

    Test-support recipe. The DFFFD flats are stacked and stored with a
    "_varAddedStack" suffix, without any stripe definition or stray light
    removal. Mostly used to test if variance is being computed correctly
    for a stack of images.

    Parameters
    ----------
    p : Primitives object
        A primitive set matching the recipe_tags.
    """
    p.prepare()
    p.checkArm()
    # p.checkND()
    p.addDQ()
    p.subtractOverscan()
    p.trimOverscan()
    p.correctImageOrientation()
    p.addVAR(read_noise=True, poisson_noise=True)
    # Creates 'DFFFD_flats' stream and leaves FDDDF flats in main stream
    p.separateFlatStreams()
    p.stackFlats(stream='DFFFD_flats', suffix='DFFFD')
    p.storeProcessedFlat(stream='DFFFD_flats', suffix='_varAddedStack')


def makeProcessedFlatDFFFF(p):
    """
    Convert raw MAROON-X flat frames into a processed DFFFF flat.

    Variant of makeProcessedFlat for datasets without FDDDF flats: the
    DDDDF flats (only fiber 5 illuminated) are kept in the main stream and
    the DFFFD flats (fibers 2, 3 and 4 illuminated) in a second stream.
    Following the legacy processing order, each stream is stacked first;
    the stacked frames are then overscan subtracted again, trimmed and
    orientation corrected. The fiber stripes of each stream are traced
    and identified, its stray light is removed, and the two streams are
    combined into a DFFFF flat (fiber 1 dark) by taking the pixel-by-pixel
    maximum. The stripe tracing is re-run on the combined frame and the
    result is stored on disk by storeProcessedFlat with a "_DFFFF_flat"
    suffix.

    Parameters
    ----------
    p : Primitives object
        A primitive set matching the recipe_tags.
    """
    p.prepare()
    p.checkArm()
    # p.checkND()
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
    p.identifyStripes(selected_fibers=[2, 3, 4, 5])
    p.defineFlatStripes(extract=True)

    p.storeProcessedFlat(suffix='_DFFFF_flat')


def makeBlaze(p):
    """
    Measure the blaze function for each fiber of a processed masterflat.

    This recipe fits a spline to the box-extracted flat spectrum of each
    fiber to model the blaze curve, then normalises each order so the peak
    equals 1. The result is stored as a BLAZE_FIBER_N extension for each
    fiber N and written to disk with a "_blaze" suffix.

    Parameters
    ----------
    p : Primitives object
        A primitive set matching the recipe_tags.
    """
    p.checkMaster()
    p.extractStripes()
    # p.boxExtraction()
    p.optimalExtraction()
    p.measureBlaze(n_knots=50)
    p.writeOutputs(suffix='_blaze')
