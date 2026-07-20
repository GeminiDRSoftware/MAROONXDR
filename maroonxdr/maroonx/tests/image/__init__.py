"""Image processing unit tests."""

import astrodata
from astropy.io import fits

import maroonx_instruments  # noqa - import is necessary for astrodata

DARK_FIBER_SETUP = ['Dark', 'Dark', 'Dark', 'Dark', 'Etalon']

# One letter per fiber type, as used in the setup code of MaroonX filenames
_FIBER_CODE = {'Dark': 'D', 'Flat lamp': 'F', 'Etalon': 'E', 'Sky': 'S', 'Target': 'O'}


def make_frame(arm, data, fiber_setup=DARK_FIBER_SETUP, nd_position=None):
    """Minimal single-extension MaroonX frame AstroData object.

    The filename has to start with a digit: a single-extension frame whose name
    starts with a letter resolves to the BUNDLE tag instead of BLUE / RED.
    """
    setup = ''.join(_FIBER_CODE[fiber] for fiber in fiber_setup)
    filename = f'00000000T000000Z_{setup}_{arm[0].lower()}_0300.fits'

    phu = fits.PrimaryHDU()
    phu.header.set('INSTRUME', 'MAROON-X')
    phu.header.set('DATALAB', 'test')
    phu.header.set('EXPTIME', 300.0)
    phu.header.set('ORIGNAME', filename)
    for number, fiber in enumerate(fiber_setup, start=1):
        phu.header.set(f'FIBER{number}', fiber)

    ext = fits.ImageHDU(data=data, name='SCI')
    ext.header.set('ARM', arm)
    ext.header.set('EXPTIME', 300.0)
    if nd_position is not None:
        ext.header.set('HIERARCH MAROONX ND POSITION', nd_position)

    ad = astrodata.create(phu, [ext])
    ad.filename = filename
    return ad
