"""
Recipes available to data with tags ['MAROONX', 'SCI'].

Default is "reduce".
"""

recipe_tags = {'MAROONX', 'SCI'}

def reduceQA(p):
    """
    Process MAROON-X science echelle spectrum with interactive QA checkpoints.

    This recipe: (1) traces and identifies the fibers and orders in a 2D
    processed flat and (2) performs both regular (aka 'box') and optimum
    extraction to produce 1D extracted spectra for 2D input spectra, with
    (3) interactive Bokeh-based spectrum visualization at key checkpoints.

    Tracing and identifying fibers and orders is done on a (preferably
    background subtracted) 2D processed flat. This step needs to be done only
    once per flat and the results can be applied to all subsequent flux
    extraction steps for other data. The routine allows to specify which fibers
    are illuminated by flat light to minimize wrong order/fiber identification.

    Box extraction is the simple summation of all spatial pixels in a given
    fiber/order combination. Optimal extraction is per default only applied to
    fibers illuminated with flat (F) and science (O) input.

    Interactive QA checkpoints use Bokeh to display spectra in a browser,
    allowing users to zoom, pan, and inspect individual orders before continuing
    reduction. The browser opens automatically at http://localhost:5006.

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
    p.addDQ()  # just placeholder until MX is in caldb
    
    p.overscanCorrect()
    p.correctImageOrientation()

    p.extractStripes()
    p.optimalExtraction(optimal_extraction_fibers=[])
    p.displaySpectra(fibers=[2, 3, 4, 5])  # QA Checkpoint 1: Verify extraction quality

    # p.getPeaksAndPolynomials(fibers=(5,), multithreading=True)
    # p.staticWavelengthSolution()
    # p.applyWavelengthSolution(fibers=(2, 3, 4), ref_fiber=5)
    # p.displaySpectra(fibers=[2, 3, 4], show_wavelength=True)  # QA Checkpoint 2: Check wavelength calibration
    # p.combineFibers()
    # p.displaySpectra(fibers=[2, 3, 4, 6], show_wavelength=True)  # QA Checkpoint 3: Inspect combined fiber spectrum
    # p.storeProcessedScience(suffix='_reducedQA')


_default = reduceQA

