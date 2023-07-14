"""
Recipes available to data with tags ['MAROONX', 'WAVECAL', 'ETALON']
Default is "reduce".

TODO: add non-default reference to static wavelength solution?
"""

recipe_tags = set(['MAROONX', 'WAVECAL', 'ETALON'])

def makeDynamicWavecal(p):
    """
    This recipe will process just etalon spectra to
    create a dynamically-by-etalon-frame updated wavelength solution reference
    for the science and sim cal fibers on the MAROON-X instrument as updated
    from a reference static wavecal in preperation for an additional drift
    calculation within a science frame during science echelle extraction.
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
    p.addVAR(read_noise=True,poisson_noise=True)
    # # get and save wavelength solution (either static reference or frame's unique sim cal solved)
    # first perform echelle extraction of fibers
    p.darkSubtraction()
    p.extractStripes()  # gets relevant flat and dark to cut out frame's spectra
    p.optimalExtraction()  # TODO: box extraction on all 5 frames
    # TODO: perform dynamic wavecal calculations on the extracted fibers
    p.storeProcessedArc(suffix='_dynamic_wavecal')
    return

_default = makeDynamicWavecal

