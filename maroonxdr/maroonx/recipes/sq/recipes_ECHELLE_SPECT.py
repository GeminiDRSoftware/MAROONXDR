"""
Recipes available to data with tags ['MAROONX', 'ECHELLE', 'SPECT']
Default is "reduce".
"""

recipe_tags = set(['MAROONX', 'ECHELLE', 'SPECT'])

def reduce(p):
    """
    This recipe processes MAROON-X echelle spectrum, up to extraction and
    stacking. (???)

    Parameters
    ----------
    p : PrimitivesCORE object
        A primitive set matching the recipe_tags.
    """

    p.prepare()
    p.addDQ()
    p.addVAR(read_noise=True)
    p.overscanCorrect()
    p.biasCorrect()
    p.ADUToElectrons()
    p.addVAR(poisson_noise=True)
    #....
    #....
    return

