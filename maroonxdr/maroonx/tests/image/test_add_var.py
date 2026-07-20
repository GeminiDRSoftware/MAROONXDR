"""Unit tests for the addVAR primitive."""

import logging

import numpy as np
import pytest

from maroonxdr.maroonx.primitives_maroonx_2D import MAROONX

from . import make_frame


# -- Fixtures ------------------------------------------------------------------
@pytest.fixture()
def data():
    """Frame data with a different photon count in every pixel."""
    return np.arange(12, dtype=np.float32).reshape(3, 4)


# -- Tests ---------------------------------------------------------------------
@pytest.mark.parametrize('arm', ['RED', 'BLUE'])
def test_addVAR(caplog, data, arm):
    """Test that addVAR adds both the read noise and the poisson contributions.

    The read noise lookup values are already quoted in variance, so they enter
    the sum as they are. The data is in DN and the gain in e-/DN, so the poisson
    contribution is the photon count converted back to a variance in DN.
    """
    caplog.set_level(logging.DEBUG)

    ad = make_frame(arm, data)
    read_noise = ad.read_noise()[0][0]
    gain = ad.gain()[0][0]

    out = MAROONX([ad]).addVAR(read_noise=True, poisson_noise=True)

    np.testing.assert_allclose(out[0][0].variance, read_noise + data / gain)
    assert len(caplog.records) > 0
    assert any('read noise variance contribution' in r.message for r in caplog.records)


@pytest.mark.parametrize('arm', ['RED', 'BLUE'])
def test_addVAR_read_noise_only(data, arm):
    """Test that the read noise alone gives a constant variance."""
    ad = make_frame(arm, data)
    read_noise = ad.read_noise()[0][0]

    out = MAROONX([ad]).addVAR(read_noise=True, poisson_noise=False)

    np.testing.assert_allclose(out[0][0].variance, read_noise)


@pytest.mark.parametrize('arm', ['RED', 'BLUE'])
def test_addVAR_poisson_noise_only(data, arm):
    """Test that the poisson noise alone gives the photon count per pixel."""
    ad = make_frame(arm, data)
    gain = ad.gain()[0][0]

    out = MAROONX([ad]).addVAR(read_noise=False, poisson_noise=True)

    np.testing.assert_allclose(out[0][0].variance, data / gain)
