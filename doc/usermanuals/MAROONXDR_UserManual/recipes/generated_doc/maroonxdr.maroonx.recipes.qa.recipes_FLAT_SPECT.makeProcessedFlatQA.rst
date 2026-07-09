makeProcessedFlatQA
===================

| **Recipe Library**: maroonxdr.maroonx.recipes.qa.recipes_FLAT_SPECT
| **Astrodata Tags**: {'MAROONX', 'CAL', 'PROCESSED', 'FLAT'}

Process MAROON-X display flat optimal extractions.

::

    Parameters
    ----------
    p : PrimitivesCORE object
        A primitive set matching the recipe_tags.

::

    def makeProcessedFlatQA(p):
        p.prepare()
        p.addDQ()  # just placeholder until MX is in caldb

        p.overscanCorrect()
        p.correctImageOrientation()

        p.extractStripes()
        p.optimalExtraction(optimal_extraction_fibers=[2, 3, 4, 5])
        p.displaySpectra(fibers=[2, 3, 4, 5])

        # measureBlaze ?
        # measureSNR ?

