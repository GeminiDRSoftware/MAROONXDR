processBundle
=============

| **Recipe Library**: maroonxdr.maroonx.recipes.sq.recipes_BUNDLE
| **Astrodata Tags**: {'MAROONX', 'BUNDLE'}

Process MAROONX observation bundles.

Red and Blue arms extensions are split before further processing.

::

    Parameters
    ----------
    p : Primitives object
        A primitive set matching the recipe_tags.

::

    def processBundle(p):
        p.splitBundle()
        p.writeOutputs()

