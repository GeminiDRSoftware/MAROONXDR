makeDarkCoefficients
====================

| **Recipe Library**: maroonxdr.maroonx.recipes.sq.recipes_DARK
| **Astrodata Tags**: {'MAROONX', 'DARK', 'CAL'}

Produce coefficient arrays z0 and z1 from pixel-by-pixel fit.

The fit uses the relationship F = z1 + z0 * log10(Texp). The input is a
list of MASTER darks created from raw dark files. The input should cover a
wide range of exposure times and can have multiple master dark files with
the same exposure time. Master darks with non-standard ND filter settings
should be avoided (they will not fit the same relationship as the other
frames).

::

    Parameters
    ----------
    p : PrimitivesCORE object
        A primitive set matching the recipe_tags.

::

    def makeDarkCoefficients(p):
        """
        Produce coefficient arrays z0 and z1 from pixel-by-pixel fit.

        The fit uses the relationship F = z1 + z0 * log10(Texp). The input is a
        list of MASTER darks created from raw dark files. The input should cover a
        wide range of exposure times and can have multiple master dark files with
        the same exposure time. Master darks with non-standard ND filter settings
        should be avoided (they will not fit the same relationship as the other
        frames).

        Parameters
        ----------
        p : PrimitivesCORE object
            A primitive set matching the recipe_tags.
        """
        p.prepare()
        p.checkArm()
        p.checkMaster()

        p.fitDarkCoefficients()
        p.storeProcessedDarkCoeff(suffix='_darkCoefficients')

