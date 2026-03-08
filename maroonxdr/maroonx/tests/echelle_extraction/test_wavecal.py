import logging
import os
from copy import deepcopy

import astrodata
import pytest

import maroonx_instruments  # noqa : import is necesary for astrodata
from maroonxdr.maroonx.primitives_maroonx_spectrum import MaroonXSpectrum


@pytest.mark.slow
@pytest.mark.preprocessed_data
@pytest.mark.parametrize('filename', ['20241124T030227Z_DEEEE_b_0030_wavecal.fits'])
def test_staticWavelengthSolution(caplog, path_to_inputs, filename):
    """
    This test checks that the static wavelength solution is correct.
    """
    caplog.set_level(logging.DEBUG)

    ad = astrodata.open(os.path.join(path_to_inputs, filename))
    p = MaroonXSpectrum([deepcopy(ad)])

    requested_fibers = (3, 4)

    adoutputs = p.staticWavelengthSolution(fibers=requested_fibers)
    out_ad = adoutputs[0]

    for fiber in [1, 2, 3, 4, 5]:
        assert hasattr(
            out_ad[0], f'WLS_STATIC_FIBER_{fiber}'
        ), f'Static solution extension for fiber {fiber} not found'

        wls = getattr(out_ad[0], f'WLS_STATIC_FIBER_{fiber}')

        if fiber in requested_fibers:
            n_orders = len(getattr(out_ad[0], f'REDUCED_ORDERS_FIBER_{fiber}'))
            assert (
                wls.shape[0] == n_orders
            ), f'Static solution for fiber {fiber} has incorrect number of orders'
        else:
            assert wls.shape == (
                1,
                1,
            ), f'Unrequested static solution for fiber {fiber} found'


@pytest.mark.slow
@pytest.mark.preprocessed_data
@pytest.mark.parametrize('filename', ['20241124T030227Z_DEEEE_b_0030_wavecal.fits'])
def test_getPeaksAndPolynomials(caplog, path_to_inputs, filename):
    """
    This test checks that PEAKS and POLY tables are set.
    """
    caplog.set_level(logging.DEBUG)

    ad = astrodata.open(os.path.join(path_to_inputs, filename))
    p = MaroonXSpectrum([deepcopy(ad)])

    p.getPeaksAndPolynomials(fibers=(4,), orders=(101,))

    out_ad = p.streams['main'][0]
    assert hasattr(out_ad[0], 'PEAKS'), 'PEAKS table should be present'
    assert hasattr(out_ad[0], 'POLY'), 'POLY table should be present'


@pytest.mark.slow
@pytest.mark.preprocessed_data
@pytest.mark.parametrize('filename', ['20241124T030227Z_DEEEE_b_0030_wavecal.fits'])
def test_fitAndApplyEtalonWls(caplog, path_to_inputs, filename):
    """
    This test checks that the dynamic wavelength solution is correct.
    """
    caplog.set_level(logging.DEBUG)

    ad = astrodata.open(os.path.join(path_to_inputs, filename))
    p = MaroonXSpectrum([deepcopy(ad)])

    requested_fibers = (3, 4)

    p.staticWavelengthSolution(fibers=requested_fibers)
    adoutputs = p.fitAndApplyEtalonWls(fibers=requested_fibers)
    out_ad = adoutputs[0]

    for fiber in [1, 2, 3, 4, 5]:
        assert hasattr(
            out_ad[0], f'WLS_DYNAMIC_FIBER_{fiber}'
        ), f'Dynamic solution extension for fiber {fiber} not found'

        wls = getattr(out_ad[0], f'WLS_DYNAMIC_FIBER_{fiber}')

        if fiber in requested_fibers:
            n_orders = len(getattr(out_ad[0], f'REDUCED_ORDERS_FIBER_{fiber}'))
            assert (
                wls.shape[0] == n_orders
            ), f'Dynamic solution for fiber {fiber} has incorrect number of orders'
        else:
            assert wls.shape == (
                1,
                1,
            ), f'Unrequested dynamic solution for fiber {fiber} found'

    # Test that updated peak data is saved with expected columns
    assert hasattr(out_ad[0], 'PEAK_DATA'), 'PEAK_DATA table should be present'
    peak_cols = set(out_ad[0].PEAK_DATA.columns)

    expected_cols = {'WAVELENGTH_BY_THAR', 'DISPERSION_MPS', 'M', 'M_FRACTION'}
    assert expected_cols.issubset(
        peak_cols
    ), f'PEAK_DATA table should contain columns: {expected_cols}'

    # Test that drift header keywords are set for requested fibers
    for fiber in requested_fibers:
        key = f'DRIFT_FIBER_{fiber}'
        assert key in out_ad[0].hdr, f'{key} header keyword should be present'
