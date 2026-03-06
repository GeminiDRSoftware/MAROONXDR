import logging
import os
from copy import deepcopy

import astrodata
import numpy as np
import pytest

import maroonx_instruments  # noqa : import is necesary for astrodata.instrument()
from maroonxdr.maroonx.primitives_maroonx_echelle import MAROONXEchelle


@pytest.mark.slow
@pytest.mark.preprocessed_data
@pytest.mark.parametrize(
    'filename',
    [
        '20241124T041907Z_SOOOE_r_0300_test_stripes.fits',
        '20241124T041907Z_SOOOE_b_0300_test_stripes.fits',
    ],
)
def test_getting_stripe_locations(caplog, path_to_inputs, filename):
    caplog.set_level(logging.DEBUG)

    ad = astrodata.open(os.path.join(path_to_inputs, filename))
    p = MAROONXEchelle([deepcopy(ad)])
    adtest = p.extractStripes()
    assert len(caplog.records) > 0
    assert any('found as associated flat' in r.message for r in caplog.records)
    # assert ad.phu.comments['REDUCTION_FLAT'] in [r.message for r in caplog.records if "found as associated flat" in r.message][0]
    # assert ad.phu.comments['REDUCTION_DARK'] in [r.message for r in caplog.records if "found as associated dark" in r.message][0]
    assert len(adtest[0][0].STRIPES.keys()) == len(ad[0].TEST_ORDERS)
    for idx, ifib in enumerate(
        sorted(adtest[0][0].STRIPES.keys(), key=lambda x: x.lower())
    ):
        np.testing.assert_allclose(
            adtest[0][0].STRIPES[ifib][str(ad[0].TEST_ORDERS[idx])].todense(),
            ad[0].STRIPES[idx],
        )
        np.testing.assert_allclose(
            adtest[0][0].F_STRIPES[ifib][str(ad[0].TEST_ORDERS[idx])].todense(),
            ad[0].F_STRIPES[idx],
        )
        np.testing.assert_allclose(
            adtest[0][0].STRIPES_MASKS[ifib][str(ad[0].TEST_ORDERS[idx])].todense(),
            ad[0].STRIPES_MASKS[idx],
        )
