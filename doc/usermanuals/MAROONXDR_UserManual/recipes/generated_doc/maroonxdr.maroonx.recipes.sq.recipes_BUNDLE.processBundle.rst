processBundle
=============

| **Recipe Library**: maroonxdr.maroonx.recipes.sq.recipes_BUNDLE
| **Astrodata Tags**: {'BUNDLE', 'MAROONX'}

Split MAROON-X observation bundles into arm-specific files.

The Red and Blue arm extensions of a raw bundle file are split into
separate AstroData objects and written to disk under their original
per-arm file names, with a header keyword added to reference the
Gemini archive bundle name. All other recipes block the BUNDLE tag,
so this recipe must run first on archive downloads.

::

    Parameters
    ----------
    p : Primitives object
        A primitive set matching the recipe_tags.

::

    def processBundle(p):
        p.splitBundle()
        p.writeOutputs()

