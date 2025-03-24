from copy import deepcopy
from pathlib import Path
import sys
import os
import logging
import pytest
import astrodata

import maroonx_instruments  # noqa : import is necesary for astrodata.instrument()
from maroonxdr.maroonx.primitives_maroonx_2D import MAROONX

# Test data should be under science_dir
science_dir = Path(__file__).parents[4] / 'science_dir'
os.chdir(science_dir)



@pytest.mark.parametrize("filename", ["20241114T181028Z_DFFFD_r_0002.fits"])
def test_nd_filter_good_series(caplog, filename):
    """
    Test that neutral density filter checker appropriately outputs the full set
    input files based on the ND filter on the sim cal fiber when all inputs are
    in agreement with the value set with the first input file.
    i.e. illumination is similar intensity as needed for good removal.
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
    assert len(adtest) == sum([test_objects[0].filter_orientation()['ND'] in
                            ad.filter_orientation()['ND'] for ad in test_objects])



@pytest.mark.parametrize("filename", ["20241114T181028Z_DFFFD_r_0002.fits"])
def test_nd_filter_subgood_series(caplog, filename):
    """
    Test that neutral density filter checker appropriately outputs the subset of
    input files based on the ND filter value set with the first input file.
    i.e. illumination is similar intensity as needed for good removal.
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
    """
       Test that neutral density filter checker appropriately IOErrors if multiple
       input files are given and the first is the only with its ND filter setting
       (i.e. somehow files were added incorrectly mid reduction).
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
    
    assert any("Only first frame found, of given, with its simcal ND filter " \
               "setting" in r.message for r in caplog.records)

if __name__ == '__main__':
    pytest.main()