"""
Recipes available to data with tags ['MAROONX', 'SCI'].

Default is "reduce".
"""

recipe_tags = {'MAROONX', 'SCI'}
blocked_tags = {'BUNDLE'}


def reduce(p):
    """
    Process MAROON-X science echelle spectrum with tracing and extraction.

    This recipe: (1) traces and identifies the fibers and orders in a 2D
    processed flat and (2) performs both regular (aka 'box') and optimum
    extraction to produce 1D extracted spectra for 2D input spectra.

    Tracing and identifying fibers and orders is done on a (preferably
    background subtracted) 2D processed flat. This step needs to be done only
    once per flat and the results can be applied to all subsequent flux
    extraction steps for other data. The routine allows to specify which fibers
    are illuminated by flat light to minimize wrong order/fiber identification.

    Box extraction is the simple summation of all spatial pixels in a given
    fiber/order combination. Optimal extraction is per default only applied to
    fibers illuminated with flat (F) and science (O) input.

    TODO: Once the Static and Dynamic wavecal recipes have been created, an
    additional set of parameters in this recipe should be added to request the
    calibration frame produced by the dynamic wavecal recipe and utilize it to
    perform a drift corrected wavelength calibration for the science frame
    fibers.

    Parameters
    ----------
    p : PrimitivesCORE object
        A primitive set matching the recipe_tags.
    """
    p.prepare()
    p.checkArm()
    p.addDQ()  # just placeholder until MX is in caldb
    p.overscanCorrect()
    p.correctImageOrientation()
    p.addVAR(read_noise=True, poisson_noise=True)
    p.extractStripes(dark_subtraction_skip_fibers=[5], straylight_removal_fibers=[5])
    p.optimalExtraction()
    p.getPeaksAndPolynomials(fibers=(5,), multithreading=True)
    p.staticWavelengthSolution()
    p.applyWavelengthSolution(fibers=(2, 3, 4), ref_fiber=5)
    p.combineFibers()
    # p.barycentricCorrection()
    p.storeProcessedScience(suffix='_reduced')


_default = reduce


def makeStripeExtractionCheck(p):
    """
    Check the stripe extraction in normal processing of a science frame.

    Parameters
    ----------
    p : PrimitivesCORE object
        A primitive set matching the recipe_tags.

    Returns
    -------
    Creates test frames with FITS-formatted stripe extractions meta-info
    (normally not saved). Unit test will independently perform stripe extraction
    and compare results.
    """
    p.prepare()
    p.checkArm()
    p.addDQ()  # just placeholder until MX is in caldb
    p.overscanCorrect()
    p.correctImageOrientation()
    p.addVAR(read_noise=True, poisson_noise=True)

    # Gets relevant flat and dark to cut out frame's spectra
    p.extractStripes(
        dark_subtraction_skip_fibers=[5],
        straylight_removal_fibers=[5],
        test_extraction=True,
    )
    p.writeOutputs(suffix='_test_stripes')


def makeSyntheticDark(p):
    """
    Construct DDDDE master darks from coefficient file for science exposures.

    Uses a linear interpolation of log(exposure time) vs. flux in empirical
    master darks to construct interpolated darks for other exposure times.

    Parameters
    ----------
    p : PrimitivesCORE object
        A primitive set matching the recipe_tags.
    """
    p.prepare()
    p.checkArm()
    p.addVAR(read_noise=True, poisson_noise=True)

    p.createSyntheticDark()
    p.storeProcessedDark(suffix='_synth_dark')


def exportReducedBundle(p):
    """
    Export reduced MAROON-X spectra from Red and Blue channels to bundle.

    Parameters
    ----------
    p : Primitives object
        A primitive set matching the recipe_tags.
    """
    p.separateArmStreams()
    p.bundleArmStreams()
    p.storeProcessedScience(suffix='_reduced')


def applyBarycentricCorrection(p):
    """
    Apply barycentric velocity correction to already-reduced MAROON-X spectra.

    Use this recipe to apply target-specific barycentric correction parameters
    after the main extraction workflow.

    Parameters
    ----------
    p : Primitives object
        A primitive set matching the recipe_tags.
    """
    p.barycentricCorrection()
    p.storeProcessedScience(suffix='_barycor')
