"""Unit tests for the fitDarkCoefficients primitive."""

import numpy as np
import pytest

from maroonxdr.maroonx.primitives_maroonx_2D import MAROONX

from . import make_frame

# A small frame keeps the per-pixel polynomial fit fast.
SHAPE = (4, 4)
# Known log-linear coefficients to recover: flux = Z1 + Z0 * log10(exptime).
Z0 = np.arange(16, dtype=np.float32).reshape(SHAPE)          # slope
Z1 = (100 - np.arange(16, dtype=np.float32)).reshape(SHAPE)  # intercept
# Powers of ten give clean log-spaced exposure times; at least five are required.
EXPOSURE_TIMES = [60., 120., 300., 600., 900., 1200., 1800.]


# -- Fixtures ------------------------------------------------------------------
@pytest.fixture()
def ad_darks():
    """Darks whose flux follows Z1 + Z0 * log10(exptime) at each exposure time."""
    ad_list = []
    for exptime in EXPOSURE_TIMES:
        flux = (Z1 + Z0 * np.log10(exptime)).astype(np.float32)
        ad = make_frame('RED', flux, nd_position=1.5)
        ad.phu['EXPTIME'] = exptime
        ad_list.append(ad)
    return ad_list


# -- Tests ---------------------------------------------------------------------
def test_fitDarkCoefficients(ad_darks):
    """The known log-linear coefficients are recovered per pixel."""
    out = MAROONX([]).fitDarkCoefficients(ad_darks)[0]

    np.testing.assert_allclose(out[0].COEFF_Z0, Z0, rtol=1e-5, atol=1e-6)
    np.testing.assert_allclose(out[0].COEFF_Z1, Z1, rtol=1e-5, atol=1e-6)
    assert out.phu.get('NCOMBINE') == len(EXPOSURE_TIMES)
    assert len(out[0].LOGEXPTIME) == len(EXPOSURE_TIMES)


def test_fitDarkCoefficients_too_few_frames(ad_darks):
    """Less than five darks cannot be fit."""
    with pytest.raises(ValueError, match='at least 5 dark frames'):
        MAROONX([]).fitDarkCoefficients(ad_darks[:4])
