"""Unit tests for the createSyntheticDarkFromCoeffs primitive."""

import numpy as np
import pytest

from maroonxdr.maroonx.primitives_maroonx_2D import MAROONX
from maroonxdr.maroonx.primitives_maroonx_echelle import MAROONXEchelle

from . import make_frame


# A small frame keeps the per-pixel polynomial fit fast.
SHAPE = (4, 4)
# Known log-linear coefficients to recover: flux = Z1 + Z0 * log10(exptime).
Z0 = np.arange(16, dtype=np.float32).reshape(SHAPE)
Z1 = (100 - np.arange(16, dtype=np.float32)).reshape(SHAPE)
# Usual dark exptimes by MX
EXPOSURE_TIMES = [60., 120., 300., 600., 900., 1200., 1800.]

# -- Fixtures ------------------------------------------------------------------
@pytest.fixture()
def dark_coeff():
    """A dark coefficients product fitted from darks with known Z0 and Z1."""
    darks = []
    for exptime in EXPOSURE_TIMES:
        flux = (Z1 + Z0 * np.log10(exptime)).astype(np.float32)
        ad = make_frame('RED', flux, nd_position=1.5)
        ad.phu['EXPTIME'] = exptime
        darks.append(ad)
    return MAROONX([]).fitDarkCoefficients(darks)[0]


# -- Tests ---------------------------------------------------------------------
def test_createSyntheticDarkFromCoeffs(dark_coeff):
    """One product per requested exposure time, named and tagged as a dark."""
    assert 'DARK_COEFF' in dark_coeff.tags
    assert 'DARK_SYNTH' not in dark_coeff.tags
    stem = dark_coeff.phu['ORIGNAME'].rsplit('_', 1)[0]

    exptimes = [120.0, 600.0]
    result = MAROONXEchelle([]).createSyntheticDarkFromCoeffs(
        [dark_coeff], exptime=exptimes)

    assert len(result) == len(exptimes)
    for ad, exptime in zip(result, exptimes):
        np.testing.assert_allclose(
            ad[0].data, Z1 + Z0 * np.log10(exptime), rtol=1e-5)

        # Exposure time field of the name rewritten, ORIGNAME kept in step.
        expected = f'{stem}_{int(exptime):04d}.fits'
        assert ad.filename == expected
        assert ad.phu['ORIGNAME'] == expected

        # EXPTIME in both headers, so descriptor and tag agree.
        assert ad.phu['EXPTIME'] == exptime
        assert ad[0].hdr['EXPTIME'] == exptime
        assert ad.exposure_time() == exptime
        assert f'{int(exptime)}s' in ad.tags

        # Coefficients gone, product is a synthetic dark and not a
        # coefficients file.
        assert not hasattr(ad[0], 'COEFF_Z0')
        assert 'DARK' in ad.tags
        assert 'DARK_SYNTH' in ad.tags
        assert 'DARK_COEFF' not in ad.tags

    # The input is left untouched.
    assert hasattr(dark_coeff[0], 'COEFF_Z0')
    assert dark_coeff.exposure_time() == EXPOSURE_TIMES[0]

    # A single exposure time may be given as a scalar.
    single = MAROONXEchelle([]).createSyntheticDarkFromCoeffs(
        [dark_coeff], exptime=120.0)
    assert len(single) == 1
    assert single[0].filename == f'{stem}_0120.fits'


def test_createSyntheticDarkFromCoeffs_requires_exptime(dark_coeff):
    """An empty or missing exptime list is rejected."""
    with pytest.raises(ValueError, match='exptime'):
        MAROONXEchelle([]).createSyntheticDarkFromCoeffs([dark_coeff])

    with pytest.raises(ValueError, match='exptime'):
        MAROONXEchelle([]).createSyntheticDarkFromCoeffs(
            [dark_coeff], exptime=[])


def test_createSyntheticDarkFromCoeffs_rejects_non_coefficients():
    """A plain dark, without the DARK_COEFF tag, is rejected."""
    dark = make_frame('RED', np.ones((4, 4), dtype=np.float32))

    with pytest.raises(ValueError, match='not a dark coefficients file'):
        MAROONXEchelle([]).createSyntheticDarkFromCoeffs(
            [dark], exptime=300)
