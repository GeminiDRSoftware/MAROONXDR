makeProcessedFlatQA
===================

| **Recipe Library**: maroonxdr.maroonx.recipes.qa.recipes_FLAT_SPECT
| **Astrodata Tags**: {'CAL', 'PROCESSED', 'MAROONX', 'FLAT'}

Inspect the extractions of a processed MAROON-X flat interactively.

QA recipe for a processed flat: the fiber and order stripes are cut out
using the flat itself as reference, optimal extraction is performed for
fibers 2 to 5, and the extracted spectra are shown at an interactive QA
checkpoint using a Bokeh server. The browser opens automatically at
localhost port 5006, allowing the user to zoom, pan and inspect
individual orders. The recipe stores no output.

::

    Parameters
    ----------
    p : Primitives object
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

