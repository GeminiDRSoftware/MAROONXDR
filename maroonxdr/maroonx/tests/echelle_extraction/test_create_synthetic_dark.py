"""Unit tests for the createSyntheticDark primitive."""

import astrodata
import numpy as np
import pytest
from astropy.io import fits
from astropy.table import Table

import maroonx_instruments  # noqa - registers the MaroonX AstroData class
from maroonxdr.maroonx.primitives_maroonx_echelle import MAROONXEchelle

from . import make_echelle_frame

# createSyntheticDark rebuilds the dark from dark = z1 + z0 * log10(exptime).
SHAPE = (4, 4)
Z0 = 2.0
Z1 = 10.0
EXPTIME = 300.0
ND_POSITION = 0.5
DARK_FIBER_SETUP = ['Dark', 'Dark', 'Dark', 'Dark', 'Etalon']


def make_dark_coeff(arm, ndfilters=(ND_POSITION,)):
    """A dark-coefficient AD with constant COEFF_Z0/Z1 and a LOGEXPTIME table."""
    phu = fits.PrimaryHDU()
    phu.header.set('INSTRUME', 'MAROON-X')

    ext = fits.ImageHDU(data=np.ones(SHAPE, dtype=np.float32), name='SCI')
    ext.header.set('ARM', arm)

    ad = astrodata.create(phu, [ext])
    ad.filename = 'dark_coeff.fits'

    ad[0].COEFF_Z0 = np.full(SHAPE, Z0, dtype=np.float32)
    ad[0].COEFF_Z1 = np.full(SHAPE, Z1, dtype=np.float32)
    ad[0].LOGEXPTIME = Table({
        'logexptime': np.log10([EXPTIME] * len(ndfilters)),
        'exptime': [EXPTIME] * len(ndfilters),
        'ndfilter': list(ndfilters),
    })
    return ad


# -- Fixtures ------------------------------------------------------------------
@pytest.fixture(params=['RED', 'BLUE'])
def ad_science(request):
    """A dark-context science frame carrying an ND filter position."""
    ad = make_echelle_frame(request.param, fiber_setup=DARK_FIBER_SETUP)
    ad[0].hdr['HIERARCH MAROONX ND POSITION'] = ND_POSITION
    return ad


# -- Tests ---------------------------------------------------------------------
def test_createSyntheticDark(ad_science):
    """The synthetic dark follows dark = z1 + z0 * log10(exptime)."""
    arm = 'BLUE' if 'BLUE' in ad_science.tags else 'RED'
    expected = Z1 + Z0 * np.log10(EXPTIME)

    dark_coeff = make_dark_coeff(arm)
    result = MAROONXEchelle([]).createSyntheticDark(
        [ad_science], dark_coeff=dark_coeff)

    assert len(result) == 1
    np.testing.assert_allclose(result[0][0].data, expected, rtol=1e-6)

    # The output shape follows the coefficient arrays.
    assert result[0][0].data.shape == SHAPE

    # Fiber keywords are stamped to a dark setup on the way out.
    for fiber in (1, 2, 3, 4):
        assert result[0].phu[f'FIBER{fiber}'] == 'Dark'
    assert result[0].phu['FIBER5'] == 'Etalon'


def test_createSyntheticDark_missing_coefficients(ad_science):
    """A dark-coefficient frame missing COEFF_Z0 is rejected."""
    arm = 'BLUE' if 'BLUE' in ad_science.tags else 'RED'

    phu = fits.PrimaryHDU()
    phu.header.set('INSTRUME', 'MAROON-X')

    ext = fits.ImageHDU(data=np.ones(SHAPE, dtype=np.float32), name='SCI')
    ext.header.set('ARM', arm)

    bad_coeff = astrodata.create(phu, [ext])
    bad_coeff.filename = 'bad_coeff.fits'

    with pytest.raises(ValueError, match='COEFF_Z0 not found'):
        MAROONXEchelle([]).createSyntheticDark([ad_science], dark_coeff=bad_coeff)
