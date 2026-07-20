"""Unit tests for the combineFlatStreams primitive."""

import logging

import numpy as np
import pytest

from maroonxdr.maroonx.primitives_maroonx_2D import MAROONX

from . import make_frame

FDDDF_FIBER_SETUP = ['Flat lamp', 'Dark', 'Dark', 'Dark', 'Flat lamp']
DFFFD_FIBER_SETUP = ['Dark', 'Flat lamp', 'Flat lamp', 'Flat lamp', 'Dark']


# -- Fixtures ------------------------------------------------------------------
def make_flat(fiber_setup, data, variance):
    """A flat frame carrying the given data and variance."""
    ad = make_frame('RED', data, fiber_setup=fiber_setup)
    ad[0].variance = variance
    return ad


@pytest.fixture()
def main_flat():
    """FDDDF flat with an ascending data pattern."""
    data = np.arange(64, dtype=np.float32).reshape(8, 8)
    return make_flat(FDDDF_FIBER_SETUP, data, data / 2)


@pytest.fixture()
def dfffd_flat():
    """DFFFD flat with a descending data pattern (complementary to main_flat)."""
    data = (63 - np.arange(64, dtype=np.float32)).reshape(8, 8)
    return make_flat(DFFFD_FIBER_SETUP, data, data / 3)


# -- Tests ---------------------------------------------------------------------
def test_combineFlatStreams(main_flat, dfffd_flat):
    """The combined flat is the per-pixel max of the two streams, data and variance."""
    expected_data = np.maximum(main_flat[0].data, dfffd_flat[0].data)
    expected_var = np.maximum(main_flat[0].variance, dfffd_flat[0].variance)

    p = MAROONX([main_flat])
    p.streams['DFFFD_flats'] = [dfffd_flat]
    out = p.combineFlatStreams(stream_2='DFFFD_flats')

    assert len(out) == 1
    np.testing.assert_array_equal(out[0].data[0], expected_data)
    np.testing.assert_allclose(out[0].variance[0], expected_var)


def test_combineFlatStreams_missing_stream(caplog, main_flat):
    """With no secondary stream the input is returned unchanged."""
    caplog.set_level(logging.DEBUG)
    original = main_flat[0].data.copy()

    p = MAROONX([main_flat])
    out = p.combineFlatStreams(stream_2='DFFFD_flats')

    assert len(out) == 1
    np.testing.assert_array_equal(out[0].data[0], original)
    assert any('does not exist so nothing to transfer' in r.message for r in caplog.records)


def test_combineFlatStreams_unexpected_lengths(caplog, main_flat, dfffd_flat):
    """A main stream that does not hold a single frame is returned unchanged."""
    caplog.set_level(logging.DEBUG)

    p = MAROONX([main_flat, main_flat])
    p.streams['DFFFD_flats'] = [dfffd_flat]
    out = p.combineFlatStreams(stream_2='DFFFD_flats')

    assert len(out) == 2
    assert any('Unexpected stream lengths' in r.message for r in caplog.records)
