makeDarkCoefficients
====================

| **Recipe Library**: maroonxdr.maroonx.recipes.sq.recipes_DARK
| **Astrodata Tags**: {'CAL', 'MAROONX', 'DARK'}

Produce coefficient arrays z0 and z1 from a pixel-by-pixel fit.

The fit is linear in log space: for each pixel, the flux in a master dark
is modeled as z1 plus z0 times log10 of the exposure time. The input is a
list of master darks created from raw dark files. The input should cover
a wide range of exposure times and can have multiple master dark files
with the same exposure time, with at least 5 master darks in total.
Master darks with non-standard ND filter settings should be avoided (they
will not fit the same relationship as the other frames). The fitted slope
and intercept arrays are stored as the COEFF_Z0 and COEFF_Z1 extensions
of a processed dark coefficients calibration.

::

    Parameters
    ----------
    p : Primitives object
        A primitive set matching the recipe_tags.

::

    def makeDarkCoefficients(p):
        p.prepare()
        p.checkArm()
        p.checkMaster()
        p.fitDarkCoefficients()
        p.storeProcessedDarkCoeff(suffix='_darkCoefficients')

