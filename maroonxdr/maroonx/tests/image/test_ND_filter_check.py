"""Tests for neutral density filter checking primitives."""

import logging
from copy import deepcopy

import astrodata
import pytest

from maroonxdr.maroonx.primitives_maroonx_2D import MAROONX


@pytest.mark.parametrize("filename", ["20241114T181028Z_DFFFD_r_0002.fits"])
def test_nd_filter_good_series(caplog, filename):
    """Test ND filter checker outputs full set when all inputs agree with first file.

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
    test_objects = [deepcopy(ad), deepcopy(ad), deepcopy(ad)]
    p = MAROONX(test_objects)
    p.prepare()
    adtest = p.checkND()

    assert len(caplog.records) > 0
    assert len(adtest) == len(test_objects)



@pytest.mark.parametrize("filename", ["20241114T181028Z_DFFFD_r_0002.fits"])
def test_nd_filter_subgood_series(caplog, filename):
    """Test ND filter checker outputs subset matching first file's ND filter value.

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
    ad_1 = deepcopy(ad)
    ad_1[0].hdr['HIERARCH MAROONX ND POSITION'] -= 12.3

    test_objects = [deepcopy(ad), deepcopy(ad), deepcopy(ad_1), deepcopy(ad_1)]
    p = MAROONX(test_objects)
    p.prepare()
    adtest = p.checkND()

    assert len(caplog.records) > 0
    assert any("Not all frames have " in r.message for r in caplog.records)
    assert len(adtest) == sum([test_objects[0].filter_orientation()['ND'] ==
                            ad.filter_orientation()['ND'] for ad in test_objects])
    assert all(test_objects[-1].filter_orientation()['ND'] !=
               ad.filter_orientation()['ND'] for ad in adtest)



@pytest.mark.parametrize("filename", ["20241114T181028Z_DFFFD_r_0002.fits"])
def test_nd_filter_bad_series(caplog, filename):
    """Test ND filter checker raises OSError when only first file has unique ND setting.

    Parameters
    ----------
    caplog : fixture
    filename : str
    """
    caplog.set_level(logging.DEBUG)

    ad = astrodata.open(filename)
    ad_1 = deepcopy(ad)
    ad_1[0].hdr['HIERARCH MAROONX ND POSITION'] -= 12.3

    test_objects = [deepcopy(ad), deepcopy(ad_1), deepcopy(ad_1), deepcopy(ad_1)]
    p = MAROONX(test_objects)
    p.prepare()

    with pytest.raises(OSError):
        p.checkND()

    assert any(
        "Only first frame found, of given, with its simcal ND filter setting"
        in r.message for r in caplog.records
    )

if __name__ == '__main__':
    pytest.main()
