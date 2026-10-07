
Add a DQ extension to the input AstroData objects.

Wrapper around the core ``addDQ`` that resolves the static bad pixel
mask for MAROON-X. With the default ``static_bpm``, the calibration
database is queried first; for every frame without a processed BPM
there, the packaged lookup for the frame's arm is used. This fallback
is permanent: MAROON-X has no archive service, so the packaged lookup
stands in for the archive-served BPM. The core primitive then flags
the bad, non-linear and saturated pixels.

Parameters
----------
adinputs : list of AstroData
    Input AstroData objects with no DQ extension.
suffix: str
    suffix to be added to output files
static_bpm: str or None
    Static bad pixel mask to apply. With "default", the calibration
    database is queried and, if it has no processed BPM for the
    frame, the packaged lookup for the frame's arm is used. With a
    file name, that file is used as given. With None, no static bad
    pixel mask is applied at all.
user_bpm: str
    Name of the bad pixel mask created by the user from flats and
    darks.  It is an optional BPM that can be added to the static one.
add_illum_mask: bool
    add illumination mask?

Returns
-------
adinputs : list of AstroData objects with a DQ extension added to them
