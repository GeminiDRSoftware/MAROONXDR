from copy import deepcopy
from pathlib import Path
import sys
import os
import logging
import pytest
import astrodata

parent_dir = Path(__file__).parents[4]
sys.path.append(str(parent_dir))
from maroonxdr.maroonx.primitives_maroonx import MAROONX
import maroonx_instruments


@pytest.mark.parametrize("DFFFD_file", ["20220725T162106Z_DFFFD_r_0001_straylight_flat.fits"])
@pytest.mark.parametrize("FDDDF_file", ["20220725T164012Z_FDDDF_r_0001_straylight_flat.fits"])
def test_stray_light_removal(caplog, DFFFD_file, FDDDF_file):
    """
    Test that removeStrayLight correctly removes the stray light across the frame.
    Requires two complimentary flat frames of a single arm that have been run
    with the makeStrayLight_check version of the FLAT_SPECT recipe
    Requires findStripes, identifyStripes, and defineFlatStripes to
    be tested. Note can create at command line with myreduce.recipename = 'makeStrayLightCheck'
    Parameters
    ----------
    caplog : fixture
    filename_r : str
    filename_b : str
    """
    caplog.set_level(logging.DEBUG)
    ad_DFFFD = astrodata.open(DFFFD_file)
    ad_FDDDF = astrodata.open(FDDDF_file)
    test_flats = [deepcopy(ad_FDDDF), deepcopy(ad_DFFFD)]
    p = MAROONX(test_flats)
    p.prepare()
    p.separateFlatStreams()
    p.findStripes()  # define stripe info to ultimately remove stray light in each stream
    p.findStripes(stream='DFFFD_flats')
    p.identifyStripes(selected_fibers='1,0,0,0,5')  # identify stripes based on MX architecture files
    p.identifyStripes(stream='DFFFD_flats', selected_fibers='0,2,3,4,0')
    p.defineFlatStripes()  # defines pixel inclusion for each flat region based on stripe ids
    p.defineFlatStripes(stream='DFFFD_flats')
    FDDDF_out = p.removeStrayLight()  # remove straylight from frame (this is why 2 partial illumination flat sets are necessary)
    DFFFD_out = p.removeStrayLight(stream='DFFFD_flats')
    assert (ad_FDDDF[0].STRAYLIGHT_DIFFERENCE == FDDDF_out[0].data - ad_FDDDF[0].data).all()
    assert (ad_DFFFD[0].STRAYLIGHT_DIFFERENCE == DFFFD_out[0].data - ad_DFFFD[0].data).all()


if __name__ == '__main__':
    pytest.main()
