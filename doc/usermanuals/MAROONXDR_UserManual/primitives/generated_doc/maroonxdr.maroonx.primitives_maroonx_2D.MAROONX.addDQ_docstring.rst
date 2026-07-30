
Add a DQ extension to the input AstroData objects.

The value of a pixel in the DQ extension will be the sum of the
following: (0=good, 1=bad pixel (found in bad pixel mask), 2=pixel is
in the non-linear regime, 4=pixel is saturated). This primitive will
trim the BPM to match the input AstroData object(s).

Parameters
----------
adinputs : list of AstroData
    Input AstroData objects with no DQ extension.
suffix: str
    suffix to be added to output files
static_bpm: str
    Name of bad pixel mask ("default" -> use default from look-up table)
    If set to None, no static_bpm will be added.
user_bpm: str
    Name of the bad pixel mask created by the user from flats and
    darks.  It is an optional BPM that can be added to the static one.
illum_mask: bool
    add illumination mask?

Returns
-------
adinputs : list of AstroData objects with a DQ extension added to them
