"""
Recipes available to data with tags ['MAROONX', 'SCI'].

Default is "reduceQA".
"""

recipe_tags = {'MAROONX', 'SCI'}
blocked_tags = {'BUNDLE'}


def reduceQA(p):
    """
    Inspect a MAROON-X science echelle spectrum at an interactive checkpoint.

    QA variant of the science reduction. The fiber and order stripes are
    cut out of the 2D frame using the traces of a processed flat and 1D
    spectra are box extracted for all fibers. The extracted spectra are
    then shown at an interactive QA checkpoint using a Bokeh server: the
    browser opens automatically at localhost port 5006, allowing the user
    to zoom, pan and inspect individual orders before continuing. The
    recipe currently stores no output; further checkpoints after
    wavelength calibration and fiber combination are present but disabled.

    Parameters
    ----------
    p : Primitives object
        A primitive set matching the recipe_tags.
    """
    p.prepare()
    p.addDQ()  # just placeholder until MX is in caldb
    p.overscanCorrect()
    p.correctImageOrientation()
    p.extractStripes()

    # empty optimal extraction fibers yield only box extraction.
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
