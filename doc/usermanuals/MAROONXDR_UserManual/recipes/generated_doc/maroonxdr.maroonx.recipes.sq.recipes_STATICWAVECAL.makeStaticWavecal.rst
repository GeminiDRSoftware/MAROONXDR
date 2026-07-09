makeStaticWavecal
=================

| **Recipe Library**: maroonxdr.maroonx.recipes.sq.recipes_STATICWAVECAL
| **Astrodata Tags**: {'MAROONX', 'ThAr', 'WAVECAL'}

Process Thorium Argon and etalon spectra for static wavelength solution.

Create a static wavelength solution reference for the science and sim cal
fibers on the MAROON-X instrument. The product of this recipe is the basis
for all wavelength calibrations on MAROON-X data (i.e. dynamical wavecals
and echelle_spect).

::

    Parameters
    ----------
    p : PrimitivesCORE object
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

