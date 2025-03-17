"""
Recipes available to data with tags ``['MAROONX', `BUNDLE`]``.
"""

recipe_tags = set(['MAROONX', 'BUNDLE', 'UNPREPARED'])


def processBundle(p):
    """
    This recipe processes MAROONX observation bundles.
    Red and Blue arms extensions are splited before further processing.

    Parameters
    ----------
    p : Primitives object
        A primitive set matching the recipe_tags.
    """
    p.splitBundle()
    p.writeOutputs()

    # for ad in p.streams['main']:
    #     ad.write(ad.phu['ORIGNAME'])


_default = processBundle
