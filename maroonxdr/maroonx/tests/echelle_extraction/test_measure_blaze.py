"""Unit tests for the measureBlaze primitive."""

import numpy as np
import pytest

from maroonxdr.maroonx.primitives_maroonx_echelle import MAROONXEchelle

from . import make_echelle_frame

# measureBlaze fits a spline to each order of BOX_REDUCED_FLAT and normalizes
# it to a per-order max of 1. A smooth Gaussian envelope is a blaze the spline
# recovers cleanly.
N_ORDERS = 3
N_PIXELS = 200
N_KNOTS = 10
FIBERS = [2, 3, 4]
ORDERS = np.array([90, 91, 92])


def make_blaze_envelope(n_orders=N_ORDERS, n_pixels=N_PIXELS):
    """Per-order Gaussian envelopes with different peak amplitudes."""
    x = np.arange(n_pixels)
    envelope = np.exp(-0.5 * ((x - n_pixels / 2) / (n_pixels / 5)) ** 2)
    # Distinct amplitude per order; normalization should erase it.
    amplitudes = np.arange(1, n_orders + 1) * 1000.0
    return (amplitudes[:, None] * envelope[None, :]).astype(np.float32)


# -- Fixtures ------------------------------------------------------------------
@pytest.fixture()
def ad_flat():
    """A processed flat carrying box-reduced Gaussian envelopes per fiber."""
    ad = make_echelle_frame('RED')
    box = make_blaze_envelope()
    for fiber in FIBERS:
        setattr(ad[0], f'BOX_REDUCED_FLAT_{fiber}', box.copy())
        setattr(ad[0], f'REDUCED_ORDERS_FIBER_{fiber}', ORDERS.copy())
    for fiber in [1, 5]:
        setattr(ad[0], f'BOX_REDUCED_FLAT_{fiber}', np.zeros((1, 1), dtype=np.float32))
    return ad


# -- Tests ---------------------------------------------------------------------
def test_measureBlaze(ad_flat):
    """Each present fiber gets a normalized blaze matching its flat's shape."""
    result = MAROONXEchelle([]).measureBlaze([ad_flat], n_knots=N_KNOTS).pop()

    for fiber in FIBERS:
        box = getattr(result[0], f'BOX_REDUCED_FLAT_{fiber}')
        blaze = getattr(result[0], f'BLAZE_FIBER_{fiber}')

        assert blaze.shape == box.shape
        assert np.all(blaze <= 1.0)
        # Every order row is normalized to a maximum of 1.
        np.testing.assert_allclose(np.nanmax(blaze, axis=1), 1.0, rtol=1e-6)


def test_measureBlaze_empty_fibers(ad_flat):
    """Fibers without a box-reduced flat get an empty placeholder blaze."""
    result = MAROONXEchelle([]).measureBlaze([ad_flat], n_knots=N_KNOTS).pop()

    for fiber in [1, 5]:
        blaze = getattr(result[0], f'BLAZE_FIBER_{fiber}')
        assert blaze.shape == (1, 1)
