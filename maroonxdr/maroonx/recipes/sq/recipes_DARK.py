"""
Recipes available to data with tags ['MAROONX', 'CAL', 'DARK'].
Default is "makeProcessedDark"
"""

recipe_tags = set(['MAROONX', 'CAL', 'DARK'])

def makeProcessedDark(p):
    """
    This recipe performs the standardization and corrections needed to convert
    the raw input dark images into a single master dark image. This output
    processed dark is stored on disk using storeProcessedDark and has a name
    equal to the name of the first input bias image with "_dark.fits" appended.

    Parameters
    ----------
    p : PrimitivesCORE object
        A primitive set matching the recipe_tags.
    """

    p.somestuff()
    return

_default = makeProcessedDark


