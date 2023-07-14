from copy import deepcopy
from pathlib import Path
import sys
import os
import logging
import pytest
import astrodata
import numpy as np

parent_dir = Path(__file__).parents[4]
sys.path.append(str(parent_dir))
from maroonxdr.maroonx.primitives_maroonx import MAROONX
import maroonx_instruments

os.chdir(parent_dir)

@pytest.mark.parametrize("filename",["science_dir/20220808T164209Z_DDDDE_r_0300.fits"])

def test_var_single(caplog, filename):
    """
    Tests the addVAR primitive works by checking if the correct arm is identified
    and if the final file contains a variance extension.
    Parameters
    ----------
    caplog : fixture
    filename : str
    Returns
    ----------
    None
    """
    caplog.set_level(logging.DEBUG)
    ad = astrodata.open(filename)
    p = MAROONX([deepcopy(ad)])
    p.prepare()
    out = p.addVAR(read_noise=True, poisson_noise=True)
    assert len(caplog.records) > 0
    assert out[0].variance is not None
    assert any("set as Red" in r.message for r in caplog.records)
