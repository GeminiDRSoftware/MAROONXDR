"""
Recipes available to data with tags ['MAROONX', 'SCI'].

Default is "reduce".
"""

recipe_tags = {'MAROONX', 'SCI'}
blocked_tags = {'BUNDLE'}


def reduce(p):
    """
    Process MAROON-X science echelle spectrum with tracing and extraction.

    This recipe: (1) traces and identifies the fibers and orders using a 2D
    processed flat, (2) performs both regular (aka 'box') and optimal
    extraction to produce 1D extracted spectra, (3) computes a drift
    corrected wavelength solution for the science fibers, and (4) combines
    the science fibers and calculates the barycentric velocity correction.
    The result is stored on disk with a "_reduced" suffix.

    Tracing and identifying fibers and orders is done on a (preferably
    background subtracted) 2D processed flat retrieved from the calibration
    database. During the stripe extraction a matching processed dark is
    subtracted from the science fibers, while the sim cal fiber gets its
    straylight removed instead.

    Box extraction is the simple summation of all spatial pixels in a given
    fiber and order combination. Optimal extraction is by default only
    applied to the science fibers 2, 3 and 4.

    The wavelength calibration fits the etalon lines of the sim cal fiber
    5, loads the static wavelength solution from a lookup file, and applies
    a drift corrected solution to the science fibers by comparison with a
    processed wavecal etalon frame retrieved from the calibration database.
    The barycentric velocity correction is computed from the exposure meter
    flux-weighted timestamps and stored in header keywords.

    Parameters
    ----------
    p : Primitives object
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
    p.getPeaksAndPolynomials(fibers=(5,))
    p.staticWavelengthSolution()
    p.applyWavelengthSolution(fibers=(2, 3, 4), ref_fiber=5)
    p.combineFibers()
    p.barycentricCorrection()
    p.storeProcessedScience(suffix='_reduced')


_default = reduce


# old recipe - set for deprecation
def makeStripeExtractionCheck(p):
    """
    Check the stripe extraction in normal processing of a science frame.

    Test-support recipe. Runs the extraction with the test option enabled so
    the stripe extractions, normally kept as sparse arrays, are saved in
    FITS-readable format (STRIPES, F_STRIPES and STRIPES_MASKS extensions).
    A unit test can then independently perform the stripe extraction and
    compare results. The output is written with a "_test_stripes" suffix.

    Parameters
    ----------
    p : Primitives object
        A primitive set matching the recipe_tags.
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
    Construct synthetic DDDDE darks for science exposures.

    The per-pixel log-linear fit stored in a processed dark coefficients
    calibration, retrieved from the calibration database, is evaluated at
    the exposure time and ND filter setting of each science frame. This
    interpolates the empirical master darks to exposure times that were not
    directly observed. The synthetic dark is stored with a "_synth_dark"
    suffix.

    Parameters
    ----------
    p : Primitives object
        A primitive set matching the recipe_tags.
    """
    p.prepare()
    p.checkArm()
    p.addVAR(read_noise=True, poisson_noise=True)
    p.createSyntheticDark()
    p.storeProcessedDark(suffix='_synth_dark')


def exportReducedBundle(p):
    """
    Bundle reduced Red and Blue arm spectra into a single output file.

    Reverses the arm split performed by processBundle: the reduced Blue and
    Red arm files of the same observation are combined into one
    multi-extension bundle, which is stored with a "_reduced" suffix.

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
    Apply barycentric velocity correction to already reduced MAROON-X spectra.

    Use this recipe to recompute the barycentric correction with target
    specific parameters (SIMBAD name, telescope coordinates, exposure meter
    zeropoints) after the main extraction workflow. The computed BERV values
    and timing information are stored as header keywords and the result is
    written with a "_barycor" suffix.

    Parameters
    ----------
    p : Primitives object
        A primitive set matching the recipe_tags.
    """
    p.barycentricCorrection()
    p.storeProcessedScience(suffix='_barycor')
