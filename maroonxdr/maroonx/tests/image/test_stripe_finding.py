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



@pytest.mark.parametrize("filename",["science_dir/20220725T164012Z_FDDDF_r_0001_FFFFF_flat.fits"])
def test_find_stripes(caplog, filename):
    """
    Test that the findStripe routine works to identify all stripes that have
    been previously found in an extracted frame as well the possible set of
    any not identified in the final extraction as being exactly those that were
    removed by identifyStripes.
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
    adtest = p.findStripes()
    all_known_stripes = np.array([])
    for keys in ad[0].STRIPES_ID.keys():
        all_known_stripes = np.append(all_known_stripes,
                                      ad[0].STRIPES_ID[keys].data)
    missing_stripes = 0
    for found_stripe in adtest[0][0].STRIPES_LOC:
        if not _stripe_in_known(all_known_stripes, found_stripe.c):
            if not (ad[0].REMOVED_STRIPES[0] == found_stripe.c).all():
                missing_stripes += 1
    assert len(caplog.records) > 0
    assert missing_stripes == 0


@pytest.mark.parametrize("filename",["science_dir/20220725T164012Z_FDDDF_r_0001_FFFFF_flat.fits"])
def test_identify_stripes(caplog, filename):
    """
    Test that the identifyStripe routine works to give order and number
    identification to stripes that have been located on MX-frames.
    Assumes that frame utilized is a previously reduced masterflat aka has 5
    illuminated fibers and an extracted STRIPES_ID extension. Checks that all
    stripes previously found in final extraction exist at identification step.
    Also check that the number of stripes removed from those inherited by
    findStripes are the same as those removed by identifyStripes
    (complementary to test in findStripes).
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
    selected_fibers = [str(ifib) for ifib in ad[0].STRIPES_FIBERS]
    p = MAROONX([deepcopy(ad)])
    p.prepare()
    p.findStripes()
    adtest = p.identifyStripes(selected_fibers=(',').join(selected_fibers))
    for ifib in selected_fibers:
        fiber = 'fiber_'+ifib
        idx = (int(ifib)-1)*6  # in the final output full order is joined
        assert fiber in adtest[0][0].STRIPES_ID.keys()
        for order in ad[0].STRIPES_ID.keys():
            if np.isnan(ad[0].STRIPES_ID[order][idx:idx+6].data).all():
                assert order not in adtest[0][0].STRIPES_ID[fiber]
            else:
                assert (adtest[0][0].STRIPES_ID[fiber][order] ==
                        ad[0].STRIPES_ID[order][idx:idx+6].data).all()
    assert len(caplog.records) > 0
    assert sum("could not be identified" in r.message for r in caplog.records) \
           == len(ad[0].REMOVED_STRIPES)



@pytest.mark.parametrize("filename", ["science_dir/20220725T164012Z_FDDDF_r_0001_FFFFF_flat.fits"])
def test_full_stripe_definition(caplog, filename):
    """
    Test that the same exact stripes are found in reference masterflat frame
    from a previous extraction
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
    selected_fibers = [str(ifib) for ifib in ad[0].STRIPES_FIBERS]
    p = MAROONX([deepcopy(ad)])
    p.prepare()
    p.findStripes()
    p.identifyStripes(selected_fibers=(',').join(selected_fibers))
    adtest = p.defineFlatStripes(extract=True)
    assert len(caplog.records) > 0
    for keys in ad[0].STRIPES_ID.keys():
        assert (ad[0].STRIPES_ID[keys] == adtest[0][0].STRIPES_ID[keys]).all()


# Local Fixtures and Helper Functions ------------------------------------------
def _stripe_in_known(A, B):
    # just a generic speedup for the array in array search
    i = 0
    j = 0
    n = len(A)
    m = len(B)
    while i < n and j < m:
        if A[i] == B[j]:
            i += 1
            j += 1
            if j == m:
                return True
        else:
            i = i - j + 1
            j = 0
    return False


if __name__ == '__main__':
    pytest.main()
