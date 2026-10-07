"""Bundle processing unit tests."""

import astrodata
import numpy as np
from astropy.io import fits

import maroonx_instruments  # noqa - import is necessary for astrodata

DARK_FIBER_SETUP = ['Dark', 'Dark', 'Dark', 'Dark', 'Etalon']

# One letter per fiber type, as used in the setup code of MaroonX filenames
_FIBER_CODE = {'Dark': 'D', 'Flat lamp': 'F', 'Etalon': 'E', 'Sky': 'S', 'Target': 'O'}


def make_arm(arm, archname, fiber_setup=DARK_FIBER_SETUP, phu_keywords=None,
             *, attach_coeffs=False):
    """Minimal single-extension MaroonX arm AstroData object.

    The filename has to start with a digit: a single-extension frame whose name
    starts with a letter resolves to the BUNDLE tag instead of BLUE / RED.

    The fiber setup selects the frame type tags (DDDDE dark by default);
    ``phu_keywords`` adds extra PHU cards (e.g. SYNTHETIC_DARK_CREATED for the
    DARK_SYNTH tag) and ``attach_coeffs`` attaches a COEFF_Z0 array to the
    extension, which is what the DARK_COEFF tag resolves from.
    """
    setup = ''.join(_FIBER_CODE[fiber] for fiber in fiber_setup)
    filename = f'00000000T000000Z_{setup}_{arm[0].lower()}_0300.fits'

    phu = fits.PrimaryHDU()
    phu.header.set('INSTRUME', 'MAROON-X')
    phu.header.set('DATALAB', 'test')
    phu.header.set('EXPTIME', 300.0)
    phu.header.set('ORIGNAME', filename)
    phu.header.set('ARCHNAME', archname)
    for number, fiber in enumerate(fiber_setup, start=1):
        phu.header.set(f'FIBER{number}', fiber)
    if phu_keywords is not None:
        for keyword, value in phu_keywords.items():
            phu.header.set(keyword, value)

    ext = fits.ImageHDU(data=np.ones((32, 32), dtype=np.float32), name='SCI')
    ext.header.set('ARM', arm)
    ext.header.set('EXPTIME', 300.0)

    ad = astrodata.create(phu, [ext])
    ad.filename = filename
    if attach_coeffs:
        ad[0].COEFF_Z0 = np.zeros((32, 32), dtype=np.float32)
    return ad
