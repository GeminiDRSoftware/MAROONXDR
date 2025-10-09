"""Tests for image orientation correction primitives."""

import logging
from copy import deepcopy

import astrodata
import numpy as np
import pytest

from maroonxdr.maroonx.primitives_maroonx_2D import MAROONX


@pytest.mark.parametrize("filename", ["20241114T181815Z_DFFFD_r_0002.fits"])
def test_correctImageOrientation_does_not_change_red_frames(caplog, filename):
    """Test that orientation does not change for raw red frames.

    Parameters
    ----------
    caplog : fixture
    filename : str

    Returns
    -------
    None
    """
    caplog.set_level(logging.DEBUG)

    ad = astrodata.open(filename)
    p = MAROONX([deepcopy(ad)])
    p.prepare()
    adtest = p.correctImageOrientation().pop()

    np.testing.assert_allclose(adtest[0].data, ad[0].data)
    assert len(caplog.records) > 0
    assert any("set as red" in r.message for r in caplog.records)


@pytest.mark.parametrize("filename", ["20241114T181959Z_DFFFD_b_0008.fits"])
def test_correctImageOrientation_flips_blue_frames(caplog, filename):
    """Test that blue frames are flipped along both axes.

    Parameters
    ----------
    caplog : fixture
    filename : str

    Returns
    -------
    None
    """
    caplog.set_level(logging.DEBUG)

    ad = astrodata.open(filename)
    p = MAROONX([deepcopy(ad)])
    p.prepare()
    adtest = p.correctImageOrientation().pop()

    np.testing.assert_allclose(adtest[0].data, np.fliplr(np.flipud(ad[0].data)))
    assert len(caplog.records) > 0
    assert any("set as blue" in r.message for r in caplog.records)


if __name__ == '__main__':
    pytest.main()
