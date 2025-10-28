"""
Recipes available to data with tags ['MAROONX', 'CAL', 'FLAT'].
"""

recipe_tags = {'MAROONX', 'CAL', 'PROCESSED', 'FLAT'}

def makeProcessedFlatQA(p):
    """
    Process MAROON-X display flat optimal extractions.

    Parameters
    ----------
    p : PrimitivesCORE object
        A primitive set matching the recipe_tags.
    """
    p.prepare()
    p.addDQ()  # just placeholder until MX is in caldb
    
    p.overscanCorrect()
    p.correctImageOrientation()

    p.extractStripes()
    p.optimalExtraction(optimal_extraction_fibers=[2, 3, 4, 5])
    p.displaySpectra(fibers=[2, 3, 4, 5])