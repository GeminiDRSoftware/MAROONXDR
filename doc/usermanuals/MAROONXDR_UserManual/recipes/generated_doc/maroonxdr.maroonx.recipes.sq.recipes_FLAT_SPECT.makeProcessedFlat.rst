makeProcessedFlat
=================

| **Recipe Library**: maroonxdr.maroonx.recipes.sq.recipes_FLAT_SPECT
| **Astrodata Tags**: {'CAL', 'MAROONX', 'FLAT'}

Convert raw MAROON-X flat frames into a single processed flat.

The input flats are separated into two streams based on their fiber
illumination pattern: FDDDF flats (fibers 1 and 5 illuminated) in the
main stream and DFFFD flats (fibers 2, 3 and 4 illuminated) in a
second stream. Each stream is stacked, its fiber stripes are traced
and identified, and its stray light is removed. The two streams are
then combined into a fully illuminated FFFFF flat by taking the
pixel-by-pixel maximum, the stripe tracing is re-run on the combined
frame, and 1D spectra are extracted for the illuminated fibers. The
result is stored on disk by storeProcessedFlat under the name of the
first input flat with "_FFFFF_flat.fits" appended.

A processed flatfield is required to perform optimal flux extraction
and to determine the blaze function. Due to the stability of the
spectrograph, one processed flatfield is typically valid for at least
a two-week period. Cross-comparisons between processed flatfields
taken months apart have not been conducted so far.

::

    Parameters
    ----------
    p : Primitives object
        A primitive set matching the recipe_tags.

::

    def makeProcessedFlat(p):
        p.prepare()
        p.checkArm()
        p.checkND()
        p.addDQ()
        p.subtractOverscan()
        p.trimOverscan()
        p.correctImageOrientation()
        p.addVAR(read_noise=True, poisson_noise=True)
        # Creates 'DFFFD_flats' stream and leaves FDDDF flats in main stream
        p.separateFlatStreams()
        p.stackFlats(suffix='FDDDF_flats')
        p.stackFlats(stream='DFFFD_flats', suffix='DFFFD')
        # Define stripe info to ultimately remove stray light in each stream
        p.findStripes()
        p.findStripes(stream='DFFFD_flats')
        # Identify stripes based on MX architecture files
        p.identifyStripes(selected_fibers=[1, 5])
        p.identifyStripes(stream='DFFFD_flats', selected_fibers=[2, 3, 4])
        # Defines pixel inclusion for each flat region based on stripe ids
        p.defineFlatStripes()
        p.defineFlatStripes(stream='DFFFD_flats')
        # Remove straylight (requires 2 partial illumination flat sets)
        p.removeStrayLight(stream='main', filter_size=19, box_size=20)
        p.removeStrayLight(stream='DFFFD_flats', filter_size=19, box_size=20)

        # Legacy patch for removeStrayLight
        # p.removeStrayLight_legacyPatch(stream='main', filter_size=19, box_size=20)
        # p.removeStrayLight_legacyPatch(stream='DFFFD_flats', filter_size=19, box_size=20)

        # Combine straylight-removed images
        p.combineFlatStreams(stream='main', stream_2='DFFFD_flats')

        # Remove second stream
        p.clearStream(stream='DFFFD_flats')
        # Re-run find/identify/define routine on combined frame
        p.findStripes()
        p.identifyStripes(selected_fibers=[1, 2, 3, 4, 5])
        p.defineFlatStripes(extract=True)

        # Perform optimal extraction on flat field to create 1D spectra
        p.extractStripes()
        p.optimalExtraction(optimal_extraction_fibers=[2, 3, 4, 5])

        p.storeProcessedFlat(suffix='_FFFFF_flat')

