import logging
import os
from copy import deepcopy

import astrodata
import numpy as np
import pytest

import maroonx_instruments  # noqa : import is necessary for astrodata.instrument()
from maroonxdr.maroonx.primitives_maroonx_echelle import MAROONXEchelle

_ALL_FIBERS = [1, 2, 3, 4, 5]


@pytest.mark.slow
@pytest.mark.preprocessed_data
@pytest.mark.parametrize(
    'filename',
    [
        '20250701T172509Z_DDDDF_b_0007_DFFFF_flat.fits',
        '20250701T171955Z_DDDDF_r_0002_DFFFF_flat.fits',
    ],
)
def test_measureBlaze_output_shape(caplog, path_to_inputs, filename):
    """
    Test that measureBlaze produces BLAZE_FIBER_{f} extensions whose shape
    matches the corresponding BOX_REDUCED_FLAT_{f} extension for every
    present fiber.
    """
    caplog.set_level(logging.DEBUG)

    ad = astrodata.open(os.path.join(path_to_inputs, filename))
    p = MAROONXEchelle([deepcopy(ad)])
    result = p.measureBlaze()
    result_ad = result[0]

    assert len(caplog.records) > 0

    for f in _ALL_FIBERS:
        box_flat = getattr(ad[0], f'BOX_REDUCED_FLAT_{f}', None)
        blaze = getattr(result_ad[0], f'BLAZE_FIBER_{f}', None)

        assert blaze is not None, f'BLAZE_FIBER_{f} extension missing from result'

        if box_flat is None or box_flat.size == 1:
            # absent fiber must produce the empty [1, 1] array
            assert blaze.shape == (
                1,
                1,
            ), f'Absent fiber {f}: expected empty shape (1, 1), got {blaze.shape}'
        else:
            # present fiber must match BOX_REDUCED_FLAT shape
            assert blaze.shape == box_flat.shape, (
                f'Fiber {f}: BLAZE_FIBER shape {blaze.shape} != '
                f'BOX_REDUCED_FLAT shape {box_flat.shape}'
            )


@pytest.mark.slow
@pytest.mark.preprocessed_data
@pytest.mark.parametrize(
    'filename',
    [
        '20250701T172509Z_DDDDF_b_0007_DFFFF_flat.fits',
        '20250701T171955Z_DDDDF_r_0002_DFFFF_flat.fits',
    ],
)
def test_measureBlaze_normalization(caplog, path_to_inputs, filename):
    """
    Test that measureBlaze produces blaze arrays that are normalized to
    have maximum value of 1.0.
    """
    caplog.set_level(logging.DEBUG)

    ad = astrodata.open(os.path.join(path_to_inputs, filename))
    p = MAROONXEchelle([deepcopy(ad)])
    result = p.measureBlaze()
    result_ad = result[0]

    for f in _ALL_FIBERS:
        blaze = getattr(result_ad[0], f'BLAZE_FIBER_{f}', None)

        assert blaze is not None, f'BLAZE_FIBER_{f} missing from result'

        assert np.all(blaze <= 1), f'BLAZE_FIBER_{f} contains values greater than 1'
