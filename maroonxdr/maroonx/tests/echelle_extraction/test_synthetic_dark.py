"""Tests for synthetic dark creation primitive."""

import numpy as np
import pytest
from astropy.io import fits
from astropy.table import Table

import astrodata
from maroonxdr.maroonx.primitives_maroonx_echelle import MAROONXEchelle


@pytest.fixture(params=['RED', 'BLUE'])
def ad_science(request):
    """Minimal science Astrodata with ND filter header."""
    
    arm = request.param
    
    phu = fits.PrimaryHDU()
    phu.header.set('INSTRUME', 'MAROON-X')
    phu.header.set('DATALAB', 'test')
    phu.header.set('EXPTIME', 300.0)
    phu.header.set('ORIGNAME', f'test_{arm[0].lower()}.fits')
    phu.header.set('FIBER1', 'Dark')
    phu.header.set('FIBER2', 'Dark')
    phu.header.set('FIBER3', 'Dark')
    phu.header.set('FIBER4', 'Dark')
    phu.header.set('FIBER5', 'Etalon')

    sci = fits.ImageHDU(data=np.ones((64, 64), dtype=np.float32), name='SCI')
    sci.header.set('ARM', arm)
    sci.header.set('EXPTIME', 300.0)
    sci.header.set('HIERARCH MAROONX ND POSITION', 0.5)

    ad = astrodata.create(phu, [sci])
    ad.filename = f'00000000T000000Z_DDDDE_{arm[0].lower()}_0300.fits'
    return ad


@pytest.fixture()
def ad_dark_coeff(ad_science):
    """Dark coefficient Astrodata with known z0, z1 values."""
    
    arm = 'RED' if 'RED' in ad_science.tags else 'BLUE'
    
    phu = fits.PrimaryHDU()
    phu.header.set('INSTRUME', 'MAROON-X')
    sci = fits.ImageHDU(data=np.ones((64, 64), dtype=np.float32), name='SCI')
    sci.header.set('ARM', arm)

    ad_coeff = astrodata.create(phu, [sci])
    ad_coeff.filename = 'dark_coeff.fits'

    ad_coeff[0].COEFF_Z0 = np.full((64, 64), 2.0, dtype=np.float32)
    ad_coeff[0].COEFF_Z1 = np.full((64, 64), 10.0, dtype=np.float32)
    ad_coeff[0].LOGEXPTIME = Table({
        'logexptime': [np.log10(300.0)],
        'exptime': [300.0],
        'ndfilter': [0.5],
    })
    return ad_coeff


def test_createSyntheticDark(ad_science, ad_dark_coeff):
    """Synthetic dark follows dark = z1 + z0 * log10(exptime)."""
    z0 = 2.0
    z1 = 10.0
    exptime = 300.0
    expected = z1 + z0 * np.log10(exptime)

    p = MAROONXEchelle([])
    result = p.createSyntheticDark([ad_science], dark_coeff=ad_dark_coeff)

    assert len(result) == 1
    np.testing.assert_allclose(result[0][0].data, expected, atol=1e-15)
