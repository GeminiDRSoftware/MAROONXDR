"""
Recipes available to data with tags ['MAROONX', 'WAVECAL', 'ECHELLE']
Default is "makeDynamicWavecal".
"""

recipe_tags = set(['MAROONX', 'WAVECAL'])

def makeDynamicWavecal(p):
    """
    This recipe processes MAROON-X 2D echelle spectra and creates the dynamic wavelength
    solution from them.  This is done in the following steps.

    1. Utilizing the relevant flat, the fibers and orders are traced and identified from the frame.
    These are extracted into sparse arrays.
    2. Using box extraction, these are converted into 1D spectra.  Box extraction is the simple summation 
    of all spatial pixels in a given fiber/order combination. The trace of the 'box' is taken from the 
    master flat field. 
    3. The extracted 1D etalon lines are fitted to determine their centroid.  This process involves 
    identifying the peaks and fitting them to a box convolved with 2 gaussians.  The locations of the 
    peak centers is stored.  The width of the peaks and the Gaussian sigmas vary very slowly along an 
    order and are modeled by low-order polynomials.
    4. The dynamic wavelength solution is computed by fitting the 1D pixel positions and the 
    wavelengths of the etalon lines using a 30 knot cubic spline.  1D pixel positions are identified by 
    comparing to the "static wavelength solution" that is stored in the caldb and is accurate up to 500 m/s.
    This allows us to calculate the "drift" of the spectrograph with time and restores accuracy to 10-20 cm/s.
    ----------
    p : PrimitivesCORE object
        A primitive set matching the recipe_tags.
    """

    p.prepare()
    p.checkArm()
    p.addDQ()  # just placeholder until MX is in caldb
    p.overscanCorrect()
    p.correctImageOrientation()
    p.addVAR(read_noise=True,poisson_noise=True)
    # # get and save wavelength solution (either static reference or frame's unique sim cal solved)
    # first perform echelle extraction of fibers
    p.extractStripes()  # gets relevant flat and dark to cut out frame's spectra
    p.boxExtraction() # extracts spectra from stripes
    p.getPeaksAndPolynomials() # fits etalon peaks and polynomials
    p.storeProcessedArc(suffix='_dynamic_wavecal')  # save reduced 1D spectra
    return

_default = makeDynamicWavecal

