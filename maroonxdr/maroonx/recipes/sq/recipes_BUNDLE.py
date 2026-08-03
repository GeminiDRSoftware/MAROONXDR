"""
Recipes available to data with tags ['MAROONX', 'BUNDLE'].

Default is "processBundle".
"""

recipe_tags = {'MAROONX', 'BUNDLE'}


def processBundle(p):
    """
    Split MAROON-X observation bundles into arm-specific files.

    The Red and Blue arm extensions of a raw bundle file are split into
    separate AstroData objects and written to disk under their original
    per-arm file names, with a header keyword added to reference the
    Gemini archive bundle name. All other recipes block the BUNDLE tag,
    so this recipe must run first on archive downloads.

    Parameters
    ----------
    p : Primitives object
        A primitive set matching the recipe_tags.
    """
    p.splitBundle()
    p.writeOutputs()


_default = processBundle





