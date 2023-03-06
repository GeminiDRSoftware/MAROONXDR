import logging
import pytest
import astrodata
from copy import deepcopy
import numpy as np
import maroonx_instruments
from maroonxdr.maroonx.primitives_maroonx import MAROONX


@pytest.mark.parametrize("filename_r", ["20200911T214124Z_DFFFD_r_0000.fits"])
@pytest.mark.parametrize("filename_b", ["20220725T162306Z_DFFFD_b_0006.fits"])
def test_checkArm_collection_and_rejection(caplog, filename_r, filename_b):
    """
    Test that first file and others of its arm-type are included
    while else are removed from set
    Parameters
    ----------
    caplog : fixture
    filename_r : str
    filename_b : str
    """
    caplog.set_level(logging.DEBUG)
    ad_red = astrodata.open(filename_r)
    ad_blue = astrodata.open(filename_b)
    test_objects = [ad_red, deepcopy(ad_red), deepcopy(ad_red),
                    ad_blue, deepcopy(ad_blue)]
    p = MAROONX(test_objects)
    p.prepare()
    out = p.checkArm()
    assert len(caplog.records) > 0
    assert any("Not all frames taken with the same camera arm"
               in r.message for r in caplog.records)
    assert all(test_objects[0].filename in ad.filename for ad in out)
    assert len(out) == sum([test_objects[0].filename in
                            ad.filename for ad in test_objects])
    assert all(test_objects[-1].filename not in ad.filename for ad in out)


@pytest.mark.parametrize("DFFFD_file", ["20220725T162306Z_DFFFD_b_0006.fits"])
@pytest.mark.parametrize("FDDDF_file", ["20220725T164341Z_FDDDF_b_0007.fits"])
def test_separating_flat_streams(caplog, DFFFD_file, FDDDF_file):
    """
    Test that seperateFlatStreams correctly separates a set of given flats by
    illuminated fibers and creates the 'DFFFD_flats' stream.
    Parameters
    ----------
    caplog : fixture
    filename_r : str
    filename_b : str
    """
    caplog.set_level(logging.DEBUG)
    ad_DFFFD = astrodata.open(DFFFD_file)
    ad_FDDDF = astrodata.open(FDDDF_file)
    test_flats = [deepcopy(ad_FDDDF), deepcopy(ad_FDDDF), deepcopy(ad_FDDDF),
                  deepcopy(ad_DFFFD), deepcopy(ad_DFFFD)]
    p = MAROONX(test_flats)
    p.prepare()
    p.separateFlatStreams()
    assert not any("Not registered as Flat" in r.message for r in caplog.records)
    assert not any("No FDDDF Flats in input list" in r.message for r in caplog.records)
    assert not any("No DFFFD Flats in input list" in r.message for r in caplog.records)
    assert len(p.streams['main']) == \
           sum([ad[0].fiber_setup() == ['Flat', 'Dark', 'Dark', 'Dark', 'Flat']
                for ad in test_flats])
    assert p.streams['DFFFD_flats'] is not None
    assert len(p.streams['DFFFD_flats']) == \
           sum([ad[0].fiber_setup() == ['Dark', 'Flat', 'Flat', 'Flat', 'Dark']
                for ad in test_flats])


@pytest.mark.parametrize("DFFFD_file", ["20220725T162306Z_DFFFD_b_0006.fits"])
@pytest.mark.parametrize("FDDDF_file", ["20220725T164341Z_FDDDF_b_0007.fits"])
def test_combining_flat_streams(caplog, DFFFD_file, FDDDF_file):
    """
    Test that combineFlatStreams correctly 'combines' a set of given two flats
    of types DFFFD and FDDDF by creating a new frame with the by-pixel max of
    the two input frames. Utilizes separateFlatStreams to create
    two-stream input.
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
    adtest = p.combineFlatStreams(stream='main', source='DFFFD_flats')
    assert not any("does not exist so nothing to transfer" in r.message for r in caplog.records)
    assert not any("Unexpected stream lengths" in r.message for r in caplog.records)
    assert len(adtest) == 1
    assert (adtest[0].data[0] == np.max([ad_DFFFD.data[0], ad_FDDDF.data[0]], axis=0)).all()


if __name__ == '__main__':
    pytest.main()
