from copy import deepcopy
from pathlib import Path
import sys
import os
import logging
import pytest
import astrodata

parent_dir = Path(__file__).parents[4]
sys.path.append(str(parent_dir))
from MAROONXDR.maroonxdr.maroonx.primitives_maroonx.primitives_maroonx_generic import MAROONX
import maroonx_instruments

os.chdir(parent_dir)

@pytest.mark.parametrize("filename",["science_dir/20220725T162106Z_DFFFD_r_0001.fits"])
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


@pytest.mark.parametrize("filename",["science_dir/20220725T162106Z_DFFFD_r_0001.fits"])
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
    ad_1.phu['HIERARCH MAROONX ND POSITION'] = \
        ad_1.phu['HIERARCH MAROONX ND POSITION'] - 10.
    test_objects = [deepcopy(ad), deepcopy(ad), deepcopy(ad_1), deepcopy(ad_1)]
    p = MAROONX(test_objects)
    p.prepare()
    adtest = p.checkND()
    assert len(caplog.records) > 0
    assert any("Not all frames have " in r.message for r in caplog.records)
    assert len(adtest) == sum([test_objects[0].filter_orientation()['ND'] in
                            ad.filter_orientation()['ND'] for ad in test_objects])
    assert all(test_objects[-1].filter_orientation()['ND'] not in
               ad.filter_orientation()['ND'] for ad in adtest)


@pytest.mark.parametrize("filename",["science_dir/20220725T162106Z_DFFFD_r_0001.fits"])
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
    ad_1.phu['HIERARCH MAROONX ND POSITION'] = \
        ad_1.phu['HIERARCH MAROONX ND POSITION'] - 10.
    test_objects = [deepcopy(ad), deepcopy(ad_1), deepcopy(ad_1), deepcopy(ad_1)]
    p = MAROONX(test_objects)
    p.prepare()
    try:
        p.checkND()
    except IOError:
        print('here')
        assert any("Only first frame found, of given, with its simcal ND filter " \
               "setting" in r.message for r in caplog.records)

if __name__ == '__main__':
    pytest.main()