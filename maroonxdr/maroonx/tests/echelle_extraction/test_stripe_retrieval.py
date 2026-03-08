import logging
import os
from copy import deepcopy
from pathlib import Path

import astrodata
import numpy as np
import pytest

import maroonx_instruments  # noqa : import is necesary for astrodata.instrument()
from maroonxdr.maroonx.primitives_maroonx_echelle import MAROONXEchelle

_CALIB_FILES = {
    'r': {
        'flat': '20241114T190714Z_DDDDF_r_0002_DFFFF_flat.fits',
        'dark': '20241124T041907Z_SOOOE_r_0300_synth_dark.fits',
    },
    'b': {
        'flat': '20241114T190714Z_DDDDF_b_0007_DFFFF_flat.fits',
        'dark': '20241124T041907Z_SOOOE_b_0300_synth_dark.fits',
    },
}


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

    arm = 'r' if '_r_' in filename else 'b'
    calib_root = Path(os.environ['DRAGONS_TEST']) / 'preprocessed_files' / 'calibrations'
    flat_path = str(calib_root / 'processed_flat' / _CALIB_FILES[arm]['flat'])
    dark_path = str(calib_root / 'processed_dark' / _CALIB_FILES[arm]['dark'])

    ad = astrodata.open(os.path.join(path_to_inputs, filename))
    p = MAROONXEchelle([deepcopy(ad)])
    adtest = p.extractStripes(
        flat=flat_path,
        dark=dark_path,
        dark_subtraction_skip_fibers=[5],
        straylight_removal_fibers=[5],
    )
    assert len(adtest[0][0].STRIPES.keys()) == len(ad[0].TEST_ORDERS)
    
    result = adtest[0][0]
    ref = ad[0]
    fibers = sorted(result.STRIPES.keys(), key=str.lower)

    for idx, fiber in enumerate(fibers):
        order = str(ref.TEST_ORDERS[idx])
        for attr in ('STRIPES', 'F_STRIPES', 'STRIPES_MASKS'):
            actual = getattr(result, attr)[fiber][order].todense()
            expected = getattr(ref, attr)[idx]
            np.testing.assert_allclose(
                actual, expected, err_msg=f'{attr} {fiber}/{order}',
            )
