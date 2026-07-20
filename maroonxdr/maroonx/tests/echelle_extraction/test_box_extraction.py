"""Unit tests for the boxExtraction primitive."""

import numpy as np
import pytest
from scipy import sparse

from maroonxdr.maroonx.primitives_maroonx_echelle import MAROONXEchelle

from . import make_echelle_frame

# Box extraction sums each column, so a full stripe of constant value gives
# value * height per column.
SIGNAL = 100.0
FLAT_SIGNAL = 1000.0
HEIGHT = 6
FIBERS = [2, 3, 4]
ORDERS = ['90', '91', '92']


def make_band(value, height=HEIGHT, ny=10, nx=8):
    """A sparse stripe: a horizontal band of constant value, ``height`` rows tall."""
    array = np.zeros((ny, nx), dtype=np.float32)
    start = (ny - height) // 2
    array[start:start + height, :] = value
    return sparse.csc_matrix(array)


def make_stripes(fibers, orders, value):
    """Build a STRIPES-style dict: {fiber_key: {order_key: sparse band}}."""
    return {
        f'fiber_{fiber}': {order: make_band(value) for order in orders}
        for fiber in fibers
    }


# -- Fixtures ------------------------------------------------------------------
@pytest.fixture()
def ad_stripes():
    """An echelle frame carrying synthetic science, flat, and mask stripes."""
    ad = make_echelle_frame('RED')
    ad[0].STRIPES = make_stripes(FIBERS, ORDERS, SIGNAL)
    ad[0].F_STRIPES = make_stripes(FIBERS, ORDERS, FLAT_SIGNAL)
    ad[0].STRIPES_MASKS = make_stripes(FIBERS, ORDERS, 1.0)
    return ad


# -- Tests ---------------------------------------------------------------------
def test_boxExtraction(ad_stripes):
    """Box-extracted flux equals the column sum of each stripe."""
    result = MAROONXEchelle([]).boxExtraction([ad_stripes]).pop()

    gain = ad_stripes.gain()[0][0]

    for fiber in FIBERS:
        box = getattr(result[0], f'BOX_REDUCED_FIBER_{fiber}')
        flat = getattr(result[0], f'BOX_REDUCED_FLAT_{fiber}')
        var = getattr(result[0], f'BOX_REDUCED_VAR_{fiber}')
        orders = getattr(result[0], f'REDUCED_ORDERS_FIBER_{fiber}')

        assert box.shape == (len(ORDERS), 8)
        np.testing.assert_allclose(box, SIGNAL * HEIGHT)
        np.testing.assert_allclose(flat, FLAT_SIGNAL * HEIGHT)
        np.testing.assert_allclose(var, SIGNAL * HEIGHT / gain)
        np.testing.assert_array_equal(orders, [float(o) for o in ORDERS])


def test_boxExtraction_absent_fibers(ad_stripes):
    """Fibers without stripes get an empty placeholder extension."""
    result = MAROONXEchelle([]).boxExtraction([ad_stripes]).pop()

    for fiber in [1, 5]:
        box = getattr(result[0], f'BOX_REDUCED_FIBER_{fiber}')
        assert box.shape == (1, 1)
