"""Echelle extraction unit tests."""

import astrodata
import numpy as np
from astropy.io import fits

import maroonx_instruments  # noqa - import is necessary for astrodata

ETALON_FIBER_SETUP = ['Dark', 'Etalon', 'Etalon', 'Etalon', 'Etalon']


def make_echelle_frame(arm='RED', fiber_setup=ETALON_FIBER_SETUP):
    """Minimal single-extension MaroonX frame for echelle extraction tests.

    Carries only what the extraction primitives read from the object itself:
    the arm (for the gain and read-noise descriptors). The stripes are attached
    separately by the test.
    """
    setup = {'Dark': 'D', 'Etalon': 'E', 'Flat lamp': 'F', 'Sky': 'S', 'Target': 'O'}
    code = ''.join(setup[fiber] for fiber in fiber_setup)
    filename = f'00000000T000000Z_{code}_{arm[0].lower()}_0300.fits'

    phu = fits.PrimaryHDU()
    phu.header.set('INSTRUME', 'MAROON-X')
    phu.header.set('DATALAB', 'test')
    phu.header.set('EXPTIME', 300.0)
    phu.header.set('ORIGNAME', filename)
    for number, fiber in enumerate(fiber_setup, start=1):
        phu.header.set(f'FIBER{number}', fiber)

    ext = fits.ImageHDU(data=np.ones((10, 8), dtype=np.float32), name='SCI')
    ext.header.set('ARM', arm)
    ext.header.set('EXPTIME', 300.0)

    ad = astrodata.create(phu, [ext])
    ad.filename = filename
    return ad
