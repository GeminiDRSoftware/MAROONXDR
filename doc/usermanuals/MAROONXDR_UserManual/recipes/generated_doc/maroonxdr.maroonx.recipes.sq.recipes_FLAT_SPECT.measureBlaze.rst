measureBlaze
============

| **Recipe Library**: maroonxdr.maroonx.recipes.sq.recipes_FLAT_SPECT
| **Astrodata Tags**: {'CAL', 'FLAT', 'MAROONX'}

Measure the blaze function for each fiber of a processed masterflat.

This primitive fits a spline to the box-extracted flat spectrum of
each fiber to model the blaze curve, then normalises each order so
the peak equals 1. The result is stored as ``BLAZE_FIBER_{f}``.

::

    Parameters
    ----------
    p : PrimitivesCORE object
        A primitive set matching the recipe_tags.

::

    def measureBlaze(p):
        """
        Measure the blaze function for each fiber of a processed masterflat.

        This primitive fits a spline to the box-extracted flat spectrum of
        each fiber to model the blaze curve, then normalises each order so
        the peak equals 1. The result is stored as ``BLAZE_FIBER_{f}``.

        Parameters
        ----------
        p : PrimitivesCORE object
            A primitive set matching the recipe_tags.
        """
        p.checkMaster()
        p.extractStripes()
        p.boxExtraction()
        # p.optimalExtraction()
        p.measureBlaze(n_knots=50)
        p.writeOutputs(suffix='_blaze')

