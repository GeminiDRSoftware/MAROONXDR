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
from MAROONXDR.maroonxdr.maroonx.primitives_maroonx.primitives_maroonx_generic import MAROONX
import maroonx_instruments

os.chdir(parent_dir)

@pytest.mark.parametrize("filename",["science_dir/20220725T162106Z_DFFFD_r_0001.fits"])
def test_correctImageOrientation_does_not_change_red_frames(caplog, filename):
    """
    Test that orientation does not change if given a raw red frame
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


@pytest.mark.parametrize("filename",["science_dir/20220725T162451Z_DFFFD_b_0006.fits"])
def test_correctImageOrientation_flips_blue_frames(caplog, filename):
    """
    Test that blue frames are flipped along both axes given a raw blue frame
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
