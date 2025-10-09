
import logging
from copy import deepcopy

import astrodata
import pytest

import maroonx_instruments  # noqa : import is necesary for astrodata.instrument()
from maroonxdr.maroonx.primitives_maroonx_echelle import MAROONXEchelle


@pytest.mark.parametrize("filename", ["20241124T062858Z_SOOOE_r_0300_reduced.fits"])
def test_getting_stripe_locations(caplog, filename):
    """
    This test checks that the stripe locations for dark subtraction are being found correctly.
    """
    caplog.set_level(logging.DEBUG)


    ad = astrodata.open(filename)
    p = MAROONXEchelle([deepcopy(ad)])
    adtest = p.attachDarkSubtraction()

    assert len(caplog.records) > 0
    assert any("Dark subtraction completed" in r.message for r in caplog.records)
    assert (adtest[0][0].DARK_SUBTRACTED == ad[0].DARK_SUBTRACTED).all()
