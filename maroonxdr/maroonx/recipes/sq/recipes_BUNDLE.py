"""Recipes available to data with tags ``['MAROONX', `BUNDLE`]``."""

recipe_tags = {'MAROONX', 'BUNDLE'}


def processBundle(p):
    """
    Process MAROONX observation bundles.

    Red and Blue arms extensions are split before further processing.

    Parameters
    ----------
    p : Primitives object
        A primitive set matching the recipe_tags.
    """
    p.splitBundle()
    p.writeOutputs()


_default = processBundle





