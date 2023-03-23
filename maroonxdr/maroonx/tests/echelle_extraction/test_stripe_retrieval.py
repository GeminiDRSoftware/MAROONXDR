import logging
import pytest
import astrodata
from copy import deepcopy
import numpy as np
import maroonx_instruments
from maroonxdr.maroonx.primitives_maroonx_echelle import MAROONXEchelle

@pytest.mark.parametrize("filename_r", ["20220808T111549Z_SOOOE_r_0300_test_stripes.fits"])
@pytest.mark.parametrize("filename_b", ["20220808T111549Z_SOOOE_b_0300_test_stripes.fits"])
def test_getting_stripe_locations(caplog, filename_r, filename_b):
    caplog.set_level(logging.DEBUG)
    for filename in [filename_b, filename_r]:
        ad = astrodata.open(filename)
        p = MAROONXEchelle([deepcopy(ad)])
        adtest = p.extractStripes()
        assert len(caplog.records) > 0
        assert any("found as associated flat" in r.message for r in caplog.records)
        assert ad.phu.comments['REDUCTION_FLAT'] in [r.message for r in caplog.records if "found as associated flat" in r.message][0]
        assert ad.phu.comments['REDUCTION_DARK'] in [r.message for r in caplog.records if "found as associated dark" in r.message][0]
        assert len(adtest[0][0].STRIPES.keys()) == len(ad[0].TEST_ORDERS)
        print('here')
        for idx,ifib in enumerate(sorted(adtest[0][0].STRIPES.keys(), key=lambda x: x.lower())):
            np.testing.assert_allclose(adtest[0][0].STRIPES[ifib][str(ad[0].TEST_ORDERS[idx])].todense(), ad[0].STRIPES[idx])
            np.testing.assert_allclose(adtest[0][0].F_STRIPES[ifib][str(ad[0].TEST_ORDERS[idx])].todense(), ad[0].F_STRIPES[idx])
            np.testing.assert_allclose(adtest[0][0].STRIPES_MASKS[ifib][str(ad[0].TEST_ORDERS[idx])].todense(), ad[0].STRIPES_MASKS[idx])