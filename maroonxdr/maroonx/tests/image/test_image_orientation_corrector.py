import logging
import pytest
import astrodata
import numpy as np
from copy import deepcopy
import maroonx_instruments
from maroonxdr.maroonx.primitives_maroonx import MAROONX


@pytest.mark.parametrize("filename",["20200911T214124Z_DFFFD_r_0000.fits"])
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


@pytest.mark.parametrize("filename",["20220725T164341Z_FDDDF_b_0007.fits"])
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
