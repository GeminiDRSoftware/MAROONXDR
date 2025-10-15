"""
Recipes available to data with tags ['MAROONX', 'WAVECAL', 'ECHELLE'].

Default is "makeDynamicWavecal".
"""

recipe_tags = {'MAROONX', 'WAVECAL'}

def makeDynamicWavecal(p):
    """
    Process MAROON-X 2D echelle spectra and create dynamic wavelength solution.

    This is done in the following steps:

    1. Utilizing the relevant flat, the fibers and orders are traced and
       identified from the frame. These are extracted into sparse arrays.
    2. Using box extraction, these are converted into 1D spectra. Box extraction
       is the simple summation of all spatial pixels in a given fiber/order
       combination. The trace of the 'box' is taken from the master flat field.
    3. The extracted 1D etalon lines are fitted to determine their centroid.
       This process involves identifying the peaks and fitting them to a box
       convolved with 2 gaussians. The locations of the peak centers is stored.
       The width of the peaks and the Gaussian sigmas vary very slowly along an
       order and are modeled by low-order polynomials.
    4. The dynamic wavelength solution is computed by fitting the 1D pixel
       positions and the wavelengths of the etalon lines using a 30 knot cubic
       spline. 1D pixel positions are identified by comparing to the "static
       wavelength solution" that is stored in the caldb and is accurate up to
       500 m/s. This allows us to calculate the "drift" of the spectrograph with
       time and restores accuracy to 10-20 cm/s.

    Parameters
    ----------
    p : PrimitivesCORE object
        A primitive set matching the recipe_tags.
    """
    p.prepare()
    p.checkArm()

    p.addDQ()  # just placeholder until MX is in caldb

    p.subtractOverscan()
    p.trimOverscan()

    p.correctImageOrientation()
    p.addVAR(read_noise=True, poisson_noise=True)
    # Get and save wavelength solution (static ref or frame's sim cal solved)
    # First perform echelle extraction of fibers
    # Gets relevant flat and dark to cut out frame's spectra
    p.extractStripes()
    # Extracts spectra from stripes
    p.boxExtraction()
    # Fits etalon peaks and polynomials
    p.getPeaksAndPolynomials(fibers=(2, 3, 4, 5), multithreading=True)

    p.staticWavelengthSolution()
    p.fitAndApplyEtalonWls()

    #p.writeOutputs(suffix='_wavecal')  # save reduced 1D spectra
    p.storeProcessedWavecal(suffix='_wavecal')  # save reduced 1D spectra

_default = makeDynamicWavecal

