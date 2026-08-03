"""Unit tests for the correctImageOrientation primitive."""

import logging

import numpy as np
import pytest

from maroonxdr.maroonx.primitives_maroonx_2D import MAROONX

from . import make_frame


# -- Fixtures ------------------------------------------------------------------
@pytest.fixture()
def data():
    """Asymmetric data, so that a flip along either axis is detectable."""
    return np.arange(12, dtype=np.float32).reshape(3, 4)


@pytest.fixture()
def mask():
    """Asymmetric mask, so that a flip along either axis is detectable."""
    return np.arange(12, dtype=np.uint16).reshape(3, 4)


# -- Tests ---------------------------------------------------------------------
def test_correctImageOrientation_no_flip_red_frames(caplog, data):
    """Test that orientation does not change for raw red frames."""
    caplog.set_level(logging.DEBUG)

    ad = make_frame('RED', data)
    p = MAROONX([ad])
    adtest = p.correctImageOrientation().pop()

    np.testing.assert_allclose(adtest[0].data, data)
    assert len(caplog.records) > 0
    assert any('set as red' in r.message for r in caplog.records)


def test_correctImageOrientation_flips_blue_frames(caplog, data):
    """Test that blue frames are flipped along both axes."""
    caplog.set_level(logging.DEBUG)

    ad = make_frame('BLUE', data)
    p = MAROONX([ad])
    adtest = p.correctImageOrientation().pop()

    np.testing.assert_allclose(adtest[0].data, np.fliplr(np.flipud(data)))
    assert len(caplog.records) > 0
    assert any('set as blue' in r.message for r in caplog.records)


def test_correctImageOrientation_flips_blue_mask(data, mask):
    """Test that the mask of blue frames is flipped along both axes."""
    ad = make_frame('BLUE', data)
    ad[0].mask = mask
    p = MAROONX([ad])
    adtest = p.correctImageOrientation().pop()

    np.testing.assert_allclose(adtest[0].mask, np.fliplr(np.flipud(mask)))


def test_correctImageOrientation_no_flip_red_mask(data, mask):
    """Test that the mask of red frames is not flipped."""
    ad = make_frame('RED', data)
    ad[0].mask = mask
    p = MAROONX([ad])
    adtest = p.correctImageOrientation().pop()

    np.testing.assert_allclose(adtest[0].mask, mask)
