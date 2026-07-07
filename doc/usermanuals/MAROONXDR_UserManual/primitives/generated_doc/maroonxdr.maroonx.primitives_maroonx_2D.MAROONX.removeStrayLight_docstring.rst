
Remove stray light from full frame images.

Remove stray light from full frame images for more accurate fiber flux
accounting. Requires the defineStripes primitive to be run prior during
recipe so INDEX_FIBER and INDEX_ORDER extensions exist to define pixel
locations across frame within fiber traces to avoid when finding stray
light.

Parameters
----------
adinputs : list of :class:`~astrodata.AstroData`
    MX flat frames. Each frame is either a DFFFD or FDDDF flat that
    has not previously had its stray light removed.

suffix : str
    Suffix to be added to output files.

box_size : int
    Pixel height and width of 'mesh_element' used in background
    identification sub-routine.

filter_size : int
    Pixel height and width of window to perform background
    identification sub-routine.

snapshot : bool
    If True, save the difference frame of removed stray light as
    the ``STRAYLIGHT_DIFFERENCE`` extension.

report : bool
    Generate PDF diagnostic report showing straylight removal stages.

Returns
-------
list of :class:`~astrodata.AstroData`
    The input frames with stray light removed from SCI, and the
    following extension added when ``snapshot=True``:

    - ``STRAYLIGHT_DIFFERENCE`` : difference frame of the removed
      stray light.
