"""Tests for frame stacking primitives."""

from copy import deepcopy

import numpy as np
import pytest
from astropy.io import fits

import astrodata
from maroonxdr.maroonx.primitives_maroonx_2D import MAROONX


@pytest.fixture(params=['RED', 'BLUE'])
def ad_flat(request):
    """Minimal MaroonX FLAT AstroData object (FDDDF fiber setup)."""
    arm = request.param
    phu = fits.PrimaryHDU()
    phu.header.set('INSTRUME', 'MAROON-X')
    phu.header.set('DATALAB', 'test')
    phu.header.set('EXPTIME', 300.0)
    phu.header.set('FIBER1', 'Flat lamp')
    phu.header.set('FIBER2', 'Dark')
    phu.header.set('FIBER3', 'Dark')
    phu.header.set('FIBER4', 'Dark')
    phu.header.set('FIBER5', 'Flat lamp')

    sci = fits.ImageHDU(data=np.ones((4400, 4400), dtype=np.float32), name='SCI')
    sci.header.set('ARM', arm)
    sci.header.set('EXPTIME', 300.0)

    ad = astrodata.create(phu, [sci])
    ad.filename = f'00000000T000000Z_FDDDF_{arm[0].lower()}_0300.fits'
    return ad


def test_stackDarks(ad_min):
    """Frames with different flux levels are scaled to the first frame before combining."""
    ad1 = ad_min
    ad2 = deepcopy(ad1)
    ad3 = deepcopy(ad1)

    levels = [100.0, 200.0, 300.0]
    for ad, level in zip([ad1, ad2, ad3], levels):
        ad[0].data = np.full_like(ad[0].data, level)
        ad.phu.set('ORIGNAME', ad.filename)

    p = MAROONX([])
    result = p.stackDarks([ad1, ad2, ad3], scale_mode='first_frame')[0]
    expected = levels[0]

    np.testing.assert_allclose(result[0].data, expected, atol=1e-15)
    assert result.phu.get('NCOMBINE') == 3


def test_stackFlats(ad_flat):
    """Frames with different flux levels are scaled to mean before combining."""
    ad1 = ad_flat
    ad2 = deepcopy(ad1)
    ad3 = deepcopy(ad1)

    levels = [100.0, 200.0, 300.0]
    for ad, level in zip([ad1, ad2, ad3], levels):
        ad[0].data = np.full_like(ad[0].data, level)
        ad.phu.set('ORIGNAME', ad.filename)

    p = MAROONX([])

    result = p.stackFlats([ad1, ad2, ad3], scale_mode='mean_frame')[0]
    expected = np.mean(levels)
    
    np.testing.assert_allclose(result[0].data, expected, atol=1e-15)
    assert result.phu.get('NCOMBINE') == 3
