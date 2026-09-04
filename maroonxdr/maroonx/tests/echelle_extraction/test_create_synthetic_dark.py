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


def make_dark_coeff(arm):
    """A dark-coefficient AD with constant COEFF_Z0/Z1 and a LOGEXPTIME table."""
    phu = fits.PrimaryHDU()
    phu.header.set('INSTRUME', 'MAROON-X')

    ext = fits.ImageHDU(data=np.ones(SHAPE, dtype=np.float32), name='SCI')
    ext.header.set('ARM', arm)

    ad = astrodata.create(phu, [ext])
    ad.filename = 'dark_coeff.fits'

    ad[0].COEFF_Z0 = np.full(SHAPE, Z0, dtype=np.float32)
    ad[0].COEFF_Z1 = np.full(SHAPE, Z1, dtype=np.float32)
    # Same columns as written by fitDarkCoefficients (no ndfilter column).
    ad[0].LOGEXPTIME = Table({
        'logexptime': [np.log10(EXPTIME)],
        'exptime': [EXPTIME],
    })
    return ad


def make_science_frame(arm, exptime, nd_position, filename):
    """A dark-context science frame with a given exposure time and ND position."""
    ad = make_echelle_frame(arm, fiber_setup=DARK_FIBER_SETUP)
    ad.phu['EXPTIME'] = exptime
    ad[0].hdr['EXPTIME'] = exptime
    ad[0].hdr['HIERARCH MAROONX ND POSITION'] = nd_position
    ad.filename = filename
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

    # Fiber keywords are stamped to a dark setup on the way out, and the
    # product tags as a synthetic dark.
    for fiber in (1, 2, 3, 4):
        assert result[0].phu[f'FIBER{fiber}'] == 'Dark'
    assert result[0].phu['FIBER5'] == 'Etalon'
    assert 'DARK_SYNTH' in result[0].tags


@pytest.mark.parametrize('arm', ['RED', 'BLUE'])
def test_createSyntheticDark_groups_same_exptime(arm):
    """Frames sharing an exposure time yield one product; ND jitter is ignored."""
    frames = [
        make_science_frame(arm, EXPTIME, ND_POSITION, 'first.fits'),
        make_science_frame(arm, EXPTIME, ND_POSITION + 0.01, 'second.fits'),
    ]
    dark_coeff = make_dark_coeff(arm)

    grouped = MAROONXEchelle([]).createSyntheticDark(
        frames, dark_coeff=dark_coeff)
    assert len(grouped) == 1
    assert grouped[0].filename == 'first.fits'

    separate = MAROONXEchelle([]).createSyntheticDark(
        frames, dark_coeff=dark_coeff, individual=True)
    assert len(separate) == 2
    assert [ad.filename for ad in separate] == ['first.fits', 'second.fits']


@pytest.mark.parametrize('arm', ['RED', 'BLUE'])
def test_createSyntheticDark_separate_exptimes(arm):
    """Frames with different exposure times each get their own product."""
    frames = [
        make_science_frame(arm, 300.0, ND_POSITION, 'short.fits'),
        make_science_frame(arm, 600.0, ND_POSITION, 'long.fits'),
    ]
    dark_coeff = make_dark_coeff(arm)

    result = MAROONXEchelle([]).createSyntheticDark(
        frames, dark_coeff=dark_coeff)

    assert [ad.filename for ad in result] == ['short.fits', 'long.fits']
    for ad, exptime in zip(result, (300.0, 600.0)):
        np.testing.assert_allclose(
            ad[0].data, Z1 + Z0 * np.log10(exptime), rtol=1e-6)


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
