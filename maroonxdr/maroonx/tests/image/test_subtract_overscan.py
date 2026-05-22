"""Tests for overscan subtraction primitive."""

import numpy as np

from maroonxdr.maroonx.primitives_maroonx_2D import MAROONX


def test_subtract_overscan(ad_min):
    """Overscan mean is subtracted from each array quadrant."""
    ad = ad_min
    signal = 100.0
    ad[0].data = np.full_like(ad[0].data, signal)

    # section where overscan is measured
    osec = ad.subtract_overscan_section()[0]
    # section where overscan is subtracted
    asec = ad.array_subtract_overscan_section()[0]

    bias = 11.0
    for sec in osec:
        ad[0].data[sec.asslice()] = signal + bias

    p = MAROONX([])
    result = p.subtractOverscan([ad])[0]
    data = result[0].data

    for osec_i in osec:
        np.testing.assert_allclose(data[osec_i.asslice()], 0.0, atol=1e-15)
