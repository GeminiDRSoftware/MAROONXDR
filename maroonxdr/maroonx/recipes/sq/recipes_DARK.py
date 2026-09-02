"""
Recipes available to data with tags ['MAROONX', 'CAL', 'DARK'].

Default is "makeProcessedDark".
"""

recipe_tags = {'MAROONX', 'CAL', 'DARK'}
blocked_tags = {'BUNDLE'}


def makeProcessedDark(p):
    """
    Convert raw MAROON-X dark frames into a single processed dark.

    The input dark frames are checked for arm and ND filter consistency,
    then stacked into one processed dark image. The result is stored on
    disk by storeProcessedDark under the name of the first input dark
    with "_dark.fits" appended.

    The background in an un-illuminated frame is very low for exposure
    times below 900s and would hardly warrant a dark subtraction on its
    own. However, most science frames are taken with the simultaneous
    calibration fiber illuminated by the FP etalon, whose extended line
    wings reach into one of the science fibers at the level of a few tens
    of counts. To remove these wings, together with any broad diffuse
    (illumination independent) background, DDDDE frames are taken in the
    daytime and stacked into a DDDDE processed dark.

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
    # legacy dark coefficients dont trim nor flip
    # p.trimOverscan()  # noqa: ERA001
    # p.correctImageOrientation()  # noqa: ERA001
    p.addVAR(read_noise=True, poisson_noise=True)
    p.stackDarks(scale_mode='first_frame', lsigma=2.0, hsigma=2.0)
    p.storeProcessedDark()


_default = makeProcessedDark


def makeDarkCoefficients(p):
    """
    Produce coefficient arrays z0 and z1 from a pixel-by-pixel fit.

    The fit is linear in log space: for each pixel, the flux in a master dark
    is modeled as z1 plus z0 times log10 of the exposure time. The input is a
    list of master darks created from raw dark files. The input should cover
    a wide range of exposure times and can have multiple master dark files
    with the same exposure time, with at least 5 master darks in total.
    Master darks with non-standard ND filter settings should be avoided (they
    will not fit the same relationship as the other frames). The fitted slope
    and intercept arrays are stored as the COEFF_Z0 and COEFF_Z1 extensions
    of a processed dark coefficients calibration.

    Parameters
    ----------
    p : Primitives object
        A primitive set matching the recipe_tags.
    """
    p.prepare()
    p.checkArm()
    p.checkMaster()
    p.fitDarkCoefficients()
    p.storeProcessedDarkCoeff(suffix='_darkCoefficients')


def makeSyntheticDarksFromCoeffs(p):
    """
    Construct synthetic DDDDE darks at given exposure times from coefficients.

    The input is a processed dark coefficients calibration (COEFF_Z0 and
    COEFF_Z1 extensions) produced by makeDarkCoefficients. For each value
    of the exptime list parameter of createSyntheticDarkFromCoeffs, the dark
    is evaluated as z1 plus z0 times log10 of the exposure time, and each
    result is stored by storeProcessedDark with a "_synth_dark" suffix.

    Parameters
    ----------
    p : Primitives object
        A primitive set matching the recipe_tags.
    """
    p.checkArm()
    p.createSyntheticDarkFromCoeffs()
    p.storeProcessedDark(suffix='_synth_dark')


# old recipe - set for deprecation
def testVARDark(p):
    """
    Produce a dark frame with an additional variance plane.

    Test-support recipe. Each input dark is processed individually (overscan
    subtracted, trimmed, orientation corrected) without stacking, and stored
    with a variance plane computed from read noise and Poisson noise. Intended
    for checking the variance computation on a single dark frame.

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
    p.storeProcessedDark(suffix='_varAdded')


# old recipe - set for deprecation
def testRegressionDark(p):
    """
    Produce a processed dark for regression comparison against legacy.

    Test-support recipe. Follows the same processing as makeProcessedDark,
    without overscan trimming nor image orientation correction to match the
    legacy dark production, and stores the result with a "_regressionDark"
    suffix for use in regression tests.

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
    # No trim overscan and no correct image orientation
    p.addVAR(read_noise=True, poisson_noise=True)
    p.stackDarks(scale_mode='first_frame', lsigma=2.0, hsigma=2.0)
    p.storeProcessedDark(suffix='_regressionDark')
