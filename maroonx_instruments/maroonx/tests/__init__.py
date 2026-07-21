"""AstroData class unit tests."""

import astrodata
import numpy as np
from astropy.io import fits

import maroonx_instruments  # noqa - import is necessary for astrodata

DARK_FIBER_SETUP = ['Dark', 'Dark', 'Dark', 'Dark', 'Etalon']
FLAT_FIBER_SETUP = ['Dark', 'Flat lamp', 'Flat lamp', 'Flat lamp', 'Dark']
ETALON_FIBER_SETUP = ['Dark', 'Etalon', 'Etalon', 'Etalon', 'Etalon']

# One letter per fiber type, as used in MaroonX filenames
_FIBER_CODE = {
    'Dark': 'D', 
    'Flat lamp': 'F', 
    'Etalon': 'E', 
    'Sky': 'S', 
    'Target': 'O'
}


def _arm_hdu(arm):
    """Minimal SCI extension for one arm."""
    ext = fits.ImageHDU(data=np.ones((32, 32), dtype=np.float32), name='SCI')
    ext.header.set('ARM', arm)
    ext.header.set('EXPTIME', 300.0)
    return ext


def _phu(fiber_setup, filename):
    """Minimal primary HDU with the keywords the adclass reads."""
    phu = fits.PrimaryHDU()
    phu.header.set('INSTRUME', 'MAROON-X')
    phu.header.set('DATALAB', 'test')
    phu.header.set('EXPTIME', 300.0)
    phu.header.set('ORIGNAME', filename)
    for number, fiber in enumerate(fiber_setup, start=1):
        phu.header.set(f'FIBER{number}', fiber)
    return phu


def make_frame(arm, fiber_setup=DARK_FIBER_SETUP):
    """Minimal single arm extension AstroData object.

    The filename has to start with a digit and follow the
    ``<timestamp>_<setup>_<arm>_<number>.fits`` convention: the arm tag and
    ``fiber_setup(short=True)`` are both derived from the filename.
    """
    setup = ''.join(_FIBER_CODE[fiber] for fiber in fiber_setup)
    filename = f'00000000T000000Z_{setup}_{arm[0].lower()}_0300.fits'

    ad = astrodata.create(
        _phu(fiber_setup, filename), 
        [_arm_hdu(arm)]
    )
    ad.filename = filename
    return ad


def make_bundle(fiber_setup=DARK_FIBER_SETUP):
    """Minimal two-extension (blue + red) MaroonX bundle AstroData object.

    The filename has to start with a letter, like the archive bundle names
    (e.g. 'N20250721M6125.fits'): a name starting with a digit is treated 
    as a split-arm frame by the arm tag logic.
    """
    filename = 'N00000000M0000.fits'

    ad = astrodata.create(
        _phu(fiber_setup, filename),
        [_arm_hdu('BLUE'), _arm_hdu('RED')]
    )
    ad.filename = filename
    return ad
