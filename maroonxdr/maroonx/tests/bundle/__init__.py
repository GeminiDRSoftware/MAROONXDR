"""Bundle processing unit tests."""

import astrodata
import numpy as np
from astropy.io import fits

import maroonx_instruments  # noqa - import is necessary for astrodata


def make_arm(arm, archname):
    """Minimal single-extension MaroonX DDDDE arm AstroData object.

    The filename has to start with a digit: a single-extension frame whose name
    starts with a letter resolves to the BUNDLE tag instead of BLUE / RED.
    """
    filename = f'00000000T000000Z_DDDDE_{arm[0].lower()}_0300.fits'

    phu = fits.PrimaryHDU()
    phu.header.set('INSTRUME', 'MAROON-X')
    phu.header.set('DATALAB', 'test')
    phu.header.set('EXPTIME', 300.0)
    phu.header.set('ORIGNAME', filename)
    phu.header.set('ARCHNAME', archname)
    for number, fiber in enumerate(['Dark', 'Dark', 'Dark', 'Dark', 'Etalon'], start=1):
        phu.header.set(f'FIBER{number}', fiber)

    ext = fits.ImageHDU(data=np.ones((32, 32), dtype=np.float32), name='SCI')
    ext.header.set('ARM', arm)
    ext.header.set('EXPTIME', 300.0)

    ad = astrodata.create(phu, [ext])
    ad.filename = filename
    return ad
