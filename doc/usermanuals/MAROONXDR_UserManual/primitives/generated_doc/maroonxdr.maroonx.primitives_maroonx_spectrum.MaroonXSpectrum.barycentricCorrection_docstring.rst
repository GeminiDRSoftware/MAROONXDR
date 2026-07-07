
Calculate barycentric velocity corrections for science observations.

This primitive computes barycentric radial velocity corrections (BERV)
for high-precision radial velocity measurements. It uses the barycorrpy
library [1]_ to calculate Earth's velocity projection toward the target
at various exposure timestamps, accounting for exposure meter
flux-weighted timing for maximum accuracy.

The primitive calculates BERV at the exposure midpoint and flux-weighted
midpoint using exposure meter data from both PC and FRD channels. It
queries SIMBAD for target coordinates or uses telescope pointing directly.
Timing corrections account for instrument-specific offsets between UTC
filename timestamps and actual exposure times.

Parameters
----------
adinputs : list of :class:`~astrodata.AstroData`
    Input AstroData objects containing 1D extracted science spectra with
    EXPOSUREMETER extension, timing information (UT_DATETIME, MJD), and
    telescope pointing data in headers.

target_name : str, optional
    Target name substring for file filtering. Only files with OBJECT
    header matching this string are processed. If None, all files are
    processed. Default is None.

simbad_target_name : str, optional
    SIMBAD-resolvable target name to override OBJECT header value. Use
    when the header target name differs from the SIMBAD catalog name.
    Only applies if target_name matches. Default is None.

use_coords : bool, optional
    If True, use telescope pointing coordinates (TELRA, TELDEC) directly
    instead of querying SIMBAD for target coordinates. Recommended when
    target is not in SIMBAD or has unreliable proper motion data.
    Default is False.

zp_pc : float, optional
    Zeropoint for PC (Precision Coupler) exposure meter channel in counts.
    If 0.0, automatically determined from median of 20 lowest readings in
    a 10-minute window around exposure. Units: counts. Default is 0.0.

zp_frd : float, optional
    Zeropoint for FRD (Fiber Refractive Diffraction) exposure meter channel
    in counts. If 0.0, automatically determined from median of 20 lowest
    readings in a 10-minute window around exposure. Units: counts.
    Default is 0.0.

start_time : str, optional
    Method for determining exposure start time. One of ``'filename'``
    (UTC from filename), ``'mjd_start'`` (telescope MJD at exposure
    start), or ``'mjd_end'`` (telescope MJD at readout end minus
    exposure time). Different methods account for varying instrument
    timing behaviors between red and blue arms. Default is ``'filename'``.

report : bool, optional
    Generate diagnostic PDF of exposure meter time series. Default is
    True.

suffix : str, optional
    Suffix to append to output filenames. Default is ``'_reduced'``.

Returns
-------
list of :class:`~astrodata.AstroData`
    Input frames with barycentric velocity corrections and exposure
    meter statistics added to the first extension header. Header
    keywords added:

    BERV values (m/s):

    - ``BERV_MIDPOINT`` : BERV at nominal exposure midpoint.
    - ``BERV_FLUXWEIGHTED_PC`` : BERV at flux-weighted midpoint (PC channel).
    - ``BERV_FLUXWEIGHTED_FRD`` : BERV at flux-weighted midpoint (FRD channel).
    - ``BERV_DIFFERENCE_PC`` : Difference between flux-weighted and nominal.
    - ``BERV_DIFFERENCE_FRD`` : Difference between flux-weighted and nominal.

    Timing information:

    - ``UTC_START`` : Corrected UTC start time (ISO format).
    - ``UTC_MIDPOINT`` : UTC midpoint time.
    - ``UTC_FLUXWEIGHTED_PC`` : Flux-weighted UTC (PC channel).
    - ``UTC_FLUXWEIGHTED_FRD`` : Flux-weighted UTC (FRD channel).
    - ``UTC_CORRECTION`` : Applied time correction in seconds.
    - ``JD_UTC_START`` : Julian date at start.
    - ``JD_UTC_MIDPOINT`` : Julian date at midpoint.
    - ``JD_UTC_FLUXWEIGHTED_PC`` : Flux-weighted JD (PC channel).
    - ``JD_UTC_FLUXWEIGHTED_FRD`` : Flux-weighted JD (FRD channel).

    Exposure meter statistics (counts):

    - ``COUNTS_PC_MIN/MAX/MEDIAN/STD`` : PC channel statistics.
    - ``COUNTS_FRD_MIN/MAX/MEDIAN/STD`` : FRD channel statistics.
    - ``COUNTS_PC_ZP`` : Applied PC zeropoint.
    - ``COUNTS_FRD_ZP`` : Applied FRD zeropoint.
    - ``SCALEFACTOR`` : Ratio of FRD to PC median counts.

    Target information:

    - ``BERV_SIMBAD_TARGET`` : Target name used for BERV calculation.

References
----------
.. [1] barycorrpy: https://github.com/shbhuk/barycorrpy
