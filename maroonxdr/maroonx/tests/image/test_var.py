"""Tests for variance extension addition primitives."""

import logging
from copy import deepcopy

import astrodata
import numpy as np
import pytest

from maroonxdr.maroonx.primitives_maroonx_2D import MAROONX


@pytest.mark.parametrize('filename', ['20241115T194624Z_DDDDE_r_0300.fits'])
def test_var_single(caplog, filename):
    """Test addVAR primitive creates variance extension correctly.

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
    out = p.addVAR(read_noise=True, poisson_noise=True)

    assert len(caplog.records) > 0
    assert out[0].variance is not None
    assert any('read noise variance contribution' in r.message for r in caplog.records)

    # check that variance array is not full of zeros
    assert not np.all(out[0].variance == 0)
