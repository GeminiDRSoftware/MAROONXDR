makeStaticWavecal
=================

| **Recipe Library**: maroonxdr.maroonx.recipes.sq.recipes_STATICWAVECAL
| **Astrodata Tags**: {'ThAr', 'MAROONX', 'WAVECAL'}

Process Thorium Argon exposures towards a static wavelength solution.

A static wavelength solution reference for the science and sim cal
fibers is the basis for all wavelength calibrations on MAROON-X data
(i.e. dynamical wavecals and science reductions). This recipe currently
performs the 2D processing and box extraction of the ThAr frames and
stores the result as a processed arc with a "_static_wavecal" suffix.
The computation of the static solution itself from the extracted lines
is not yet implemented; the pipeline relies on the static solution
distributed as a lookup file.

::

    Parameters
    ----------
    p : Primitives object
        A primitive set matching the recipe_tags.

::

    def makeStaticWavecal(p):
        p.prepare()
        p.checkArm()
        p.addDQ()  # just placeholder until MX is in caldb
        p.overscanCorrect()
        p.correctImageOrientation()
        p.addVAR(read_noise=True, poisson_noise=True)
        # Get and save wavelength solution (static ref or frame's sim cal solved)
        # First perform echelle extraction of fibers
        # Gets relevant flat and dark to cut out frame's spectra
        p.extractStripes()
        p.boxExtraction()
        # TODO: second perform static wavecal calculations on the extracted fibers
        #
        p.storeProcessedArc(suffix='_static_wavecal')

