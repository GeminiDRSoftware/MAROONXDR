import logging
import pytest
import astrodata
from copy import deepcopy
import numpy as np
import maroonx_instruments
from maroonxdr.maroonx.primitives_maroonx_echelle import MAROONXEchelle


@pytest.mark.parametrize("filename_r", ["20220808T111549Z_SOOOE_r_0300_reduced.fits"])
@pytest.mark.parametrize("filename_b", ["20220808T111549Z_SOOOE_b_0300_reduced.fits"])
def test_optimal_extracting_science_data(caplog, filename_r, filename_b):
    """
    Test that new optimal extraction and box extraction is equal to
    a previously reduced extraction of the same intial dataset
    for a standard science frame  (aka 3 science fibers, sim cal and sky fiber).
    As is standard for the echelle extraction, all extensions should exist even
    if not called to be populated by data
    (e.g. optimal extraction of simcal fiber).
    Parameters
    ----------
    caplog : fixture
    filename_r : str
    filename_b : str
    """
    caplog.set_level(logging.DEBUG)
    for filename in [filename_b, filename_r]:
        ad = astrodata.open(filename)
        p = MAROONXEchelle([deepcopy(ad)])
        p.extractStripes()
        adtest = p.optimalExtraction()
        assert len(caplog.records) > 0
        assert any("extracted" in r.message for r in caplog.records)
        # first check that all orders were found for all fibers
        assert (ad[0].REDUCED_ORDERS_FIBER_1 == adtest[0][0].REDUCED_ORDERS_FIBER_1).all()
        assert (ad[0].REDUCED_ORDERS_FIBER_2 == adtest[0][0].REDUCED_ORDERS_FIBER_2).all()
        assert (ad[0].REDUCED_ORDERS_FIBER_3 == adtest[0][0].REDUCED_ORDERS_FIBER_3).all()
        assert (ad[0].REDUCED_ORDERS_FIBER_4 == adtest[0][0].REDUCED_ORDERS_FIBER_4).all()
        assert (ad[0].REDUCED_ORDERS_FIBER_5 == adtest[0][0].REDUCED_ORDERS_FIBER_5).all()
        # then check optimal data sets
        assert np.array_equal(ad[0].OPTIMAL_REDUCED_FIBER_1, adtest[0][0].OPTIMAL_REDUCED_FIBER_1, equal_nan=True)
        assert np.array_equal(ad[0].OPTIMAL_REDUCED_ERR_1, adtest[0][0].OPTIMAL_REDUCED_ERR_1, equal_nan=True)
        assert np.array_equal(ad[0].OPTIMAL_REDUCED_FIBER_2, adtest[0][0].OPTIMAL_REDUCED_FIBER_2, equal_nan=True)
        assert np.array_equal(ad[0].OPTIMAL_REDUCED_ERR_2, adtest[0][0].OPTIMAL_REDUCED_ERR_2, equal_nan=True)
        assert np.array_equal(ad[0].OPTIMAL_REDUCED_FIBER_3, adtest[0][0].OPTIMAL_REDUCED_FIBER_3, equal_nan=True)
        assert np.array_equal(ad[0].OPTIMAL_REDUCED_ERR_3, adtest[0][0].OPTIMAL_REDUCED_ERR_3, equal_nan=True)
        assert np.array_equal(ad[0].OPTIMAL_REDUCED_FIBER_4, adtest[0][0].OPTIMAL_REDUCED_FIBER_4, equal_nan=True)
        assert np.array_equal(ad[0].OPTIMAL_REDUCED_ERR_4, adtest[0][0].OPTIMAL_REDUCED_ERR_4, equal_nan=True)
        assert np.array_equal(ad[0].OPTIMAL_REDUCED_FIBER_5, adtest[0][0].OPTIMAL_REDUCED_FIBER_5, equal_nan=True)
        assert np.array_equal(ad[0].OPTIMAL_REDUCED_ERR_5, adtest[0][0].OPTIMAL_REDUCED_ERR_5, equal_nan=True)
        # then check box reductions
        assert np.array_equal(ad[0].BOX_REDUCED_FIBER_1, adtest[0][0].BOX_REDUCED_FIBER_1, equal_nan=True)
        assert np.array_equal(ad[0].BOX_REDUCED_ERR_1, adtest[0][0].BOX_REDUCED_ERR_1, equal_nan=True)
        assert np.array_equal(ad[0].BOX_REDUCED_FIBER_2, adtest[0][0].BOX_REDUCED_FIBER_2, equal_nan=True)
        assert np.array_equal(ad[0].BOX_REDUCED_ERR_2, adtest[0][0].BOX_REDUCED_ERR_2, equal_nan=True)
        assert np.array_equal(ad[0].BOX_REDUCED_FIBER_3, adtest[0][0].BOX_REDUCED_FIBER_3, equal_nan=True)
        assert np.array_equal(ad[0].BOX_REDUCED_ERR_3, adtest[0][0].BOX_REDUCED_ERR_3, equal_nan=True)
        assert np.array_equal(ad[0].BOX_REDUCED_FIBER_4, adtest[0][0].BOX_REDUCED_FIBER_4, equal_nan=True)
        assert np.array_equal(ad[0].BOX_REDUCED_ERR_4, adtest[0][0].BOX_REDUCED_ERR_4, equal_nan=True)
        assert np.array_equal(ad[0].BOX_REDUCED_FIBER_5, adtest[0][0].BOX_REDUCED_FIBER_5, equal_nan=True)
        assert np.array_equal(ad[0].BOX_REDUCED_ERR_5, adtest[0][0].BOX_REDUCED_ERR_5, equal_nan=True)
        # finally check the bad pixel masks for the extracted spectra
        assert (ad[0].BPM_FIBER_1 == adtest[0][0].BPM_FIBER_1).all()
        assert (ad[0].BPM_FIBER_2 == adtest[0][0].BPM_FIBER_2).all()
        assert (ad[0].BPM_FIBER_3 == adtest[0][0].BPM_FIBER_3).all()
        assert (ad[0].BPM_FIBER_4 == adtest[0][0].BPM_FIBER_4).all()
        assert (ad[0].BPM_FIBER_5 == adtest[0][0].BPM_FIBER_5).all()


if __name__ == '__main__':
    pytest.main()
