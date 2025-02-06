"""
Recipes available to data with tags ['MAROONX', 'WAVECAL', 'ThAr']
Default is "reduce".

"""

recipe_tags = set(['MAROONX', 'WAVECAL', 'ThAr'])

def makeStaticWavecal(p):
    """
    This recipe will process Thorium Argon spectra and etalon spectra to
    create a static wavelength solution reference for the science and sim cal
    fibers on the MAROON-X instrument. The product of this recipe is the bases
    for all wavelength calibrations on MAROON-X data (i.e. dynamical wavecals
    and echelle_spect).
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
    # get and save wavelength solution (either static reference or frame's unique sim cal solved)
    # first perform echelle extraction of fibers
    p.extractStripes()  # gets relevant flat and dark to cut out frame's spectra
    p.boxExtraction()
    # TODO: second perform static wavecal calculations on the extracted fibers
    # 
    p.storeProcessedArc(suffix='_static_wavecal')
    return

_default = makeStaticWavecal

