import logging
import os
from copy import deepcopy
from pathlib import Path

import astrodata
import numpy as np
import pytest

import maroonx_instruments  # noqa : import is necesary for astrodata.instrument()
from maroonxdr.maroonx.primitives_maroonx_2D import MAROONX

# Test data should be under science_dir
science_dir = Path(__file__).parents[4] / 'science_dir'
os.chdir(science_dir)


@pytest.mark.parametrize('filename', ['20241115T194624Z_DDDDE_r_0300.fits'])
def test_var_single(caplog, filename):
    """
    Tests the addVAR primitive works by checking if the correct arm is identified
    and if the final file contains a variance extension.

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
