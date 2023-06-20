import logging
import pytest
import astrodata
from copy import deepcopy
import numpy as np
import maroonx_instruments
from maroonxdr.maroonx.primitives_maroonx_echelle import MAROONXEchelle

@pytest.mark.parametrize("filename", ["./maroonxdr/maroonx/tests/echelle_extraction/20220808T111549Z_SOOOE_r_0300_reduced.fits"])
def test_getting_stripe_locations(caplog, filename):
    """
    This test checks that the stripe locations for dark subtraction are being found correctly.
    """
    caplog.set_level(logging.DEBUG)
    ad = astrodata.open(filename)
    p = MAROONXEchelle([deepcopy(ad)])
    adtest = p.darkSubtraction()
    print('here')
    assert len(caplog.records) > 0
    assert any("found as associated " in r.message for r in caplog.records)
    assert (adtest[0][0].DARK_SUBTRACTED == ad[0].DARK_SUBTRACTED).all()