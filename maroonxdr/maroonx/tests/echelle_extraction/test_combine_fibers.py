"""Unit tests for the combineFibers primitive.

Usage:
    python -m maroonxdr.maroonx.tests.echelle_extraction.test_combine_fibers --create-inputs
"""

import os
from copy import deepcopy

import numpy as np
import pytest
from gempy.utils import logutils

import astrodata
import maroonx_instruments  # noqa - registers the MaroonX AstroData class
from maroonxdr.maroonx.primitives_maroonx_spectrum import MaroonXSpectrum

from . import make_echelle_frame

# Synthetic spectra: constant flux and variance per fiber. Fibers 2 and 4 are
# scaled to the flux level of fiber 3 (the reference) before combination, so
# the combined spectrum must equal the fiber 3 level everywhere.
NPIX = 256
ORDERS = [90.0, 91.0, 92.0]
FLUX = {2: 200.0, 3: 100.0, 4: 50.0}
VAR = {2: 4.0, 3: 2.0, 4: 1.0}
# Scaling divides each fiber's variance by its flux ratio to fiber 3, so all
# three fibers carry variance 2 and weight 1/2 into the combination.
COMBINED_VAR = 1.0 / 1.5

# Blessed reduced science products, per arm.
datasets = [
    '20250717T144308Z_SOOOE_b_0300_reduced.fits',
    '20250717T144308Z_SOOOE_r_0300_reduced.fits',
]

# Parameters of the blessed run (preprocess/science.py); the rest are defaults.
COMBINE_PARAMS = {'max_clips': 20, 'report': False}


def make_wls():
    """A strictly increasing wavelength grid per order, shape (orders, NPIX)."""
    rows = []
    for i in range(len(ORDERS)):
        start = 500.0 + 10.0 * i
        rows.append(np.linspace(start, start + 5.0, NPIX))
    return np.vstack(rows)


def make_ad_fibers(arm, flux=FLUX, var=VAR):
    """An echelle frame carrying extracted spectra for fibers 2, 3, and 4."""
    ad = make_echelle_frame(arm)
    for fiber in (2, 3, 4):
        setattr(ad[0], f'OPTIMAL_REDUCED_FIBER_{fiber}',
                np.full((len(ORDERS), NPIX), flux[fiber]))
        setattr(ad[0], f'OPTIMAL_REDUCED_VAR_{fiber}',
                np.full((len(ORDERS), NPIX), var[fiber]))
        setattr(ad[0], f'WLS_SIMULTANEOUS_FIBER_{fiber}', make_wls())
        setattr(ad[0], f'REDUCED_ORDERS_FIBER_{fiber}', ORDERS.copy())
    return ad


# -- Tests ---------------------------------------------------------------------
@pytest.mark.parametrize('arm', ['RED', 'BLUE'])
def test_combineFibers(arm):
    """Scaled constant fibers combine to the reference-fiber flux level."""
    ad = make_ad_fibers(arm)
    result = MaroonXSpectrum([]).combineFibers([ad], report=False).pop()

    flux = result[0].OPTIMAL_REDUCED_FIBER_6
    var = result[0].OPTIMAL_REDUCED_VAR_6

    assert flux.shape == (len(ORDERS), NPIX)
    # The known-bad pixel (196 in the blue arm) is masked and re-interpolated,
    # so the flux is flat at the reference level everywhere in both arms.
    np.testing.assert_allclose(flux, FLUX[3], rtol=1e-9)

    expected_var = np.full((len(ORDERS), NPIX), COMBINED_VAR)
    if arm == 'BLUE':
        expected_var[:, 196] = 1e6  # masked pixel carries no information
    np.testing.assert_allclose(var, expected_var, rtol=1e-12)

    # Wavelength grid and orders are taken from the reference fiber
    np.testing.assert_array_equal(
        result[0].WLS_SIMULTANEOUS_FIBER_6,
        result[0].WLS_SIMULTANEOUS_FIBER_3,
    )
    np.testing.assert_array_equal(result[0].REDUCED_ORDERS_FIBER_6, ORDERS)


def test_combineFibers_rejects_outlier():
    """A single-pixel outlier is kappa-sigma clipped from the combination."""
    spike = 50
    flux = {2: 100.0, 3: 100.0, 4: 100.0}
    var = {2: 1.0, 3: 1.0, 4: 1.0}

    def ad_with_spike():
        ad = make_ad_fibers('RED', flux=flux, var=var)
        ad[0].OPTIMAL_REDUCED_FIBER_4[0, spike] = 300.0
        return ad

    # At the default threshold the outlier is rejected: the combination at the
    # spike pixel falls back to fibers 2 and 3.
    clipped = MaroonXSpectrum([]).combineFibers([ad_with_spike()], report=False).pop()
    flux6 = clipped[0].OPTIMAL_REDUCED_FIBER_6
    np.testing.assert_allclose(flux6[0, spike], flux[3], rtol=1e-9)

    # With the threshold effectively disabled the outlier leaks through.
    unclipped = MaroonXSpectrum([]).combineFibers(
        [ad_with_spike()], kappa_sigma=1000.0, report=False
    ).pop()
    flux6 = unclipped[0].OPTIMAL_REDUCED_FIBER_6
    assert flux6[0, spike] > 150.0


@pytest.mark.preprocessed_data
@pytest.mark.regression
@pytest.mark.parametrize('filename', datasets)
def test_combineFibers_real_science(filename, path_to_inputs, change_working_dir):
    """Fiber combination of a real science frame is stable against the ref."""
    ad = astrodata.open(os.path.join(path_to_inputs, filename))

    # Stored output from the reference run
    ref = {
        ext: getattr(ad[0], ext)
        for ext in ('OPTIMAL_REDUCED_FIBER_6', 'OPTIMAL_REDUCED_VAR_6',
                    'WLS_SIMULTANEOUS_FIBER_6', 'REDUCED_ORDERS_FIBER_6')
    }

    with change_working_dir():
        logutils.config(file_name=f'log_{filename.replace(".fits", "")}.txt')
        ad_out = MaroonXSpectrum([]).combineFibers(
            [deepcopy(ad)], max_clips=20, report=False
        )[0]

    for ext in ('OPTIMAL_REDUCED_FIBER_6', 'OPTIMAL_REDUCED_VAR_6'):
        np.testing.assert_allclose(
            getattr(ad_out[0], ext), ref[ext],
            rtol=1e-6, atol=1e-8, err_msg=ext,
        )
    np.testing.assert_array_equal(
        ad_out[0].WLS_SIMULTANEOUS_FIBER_6, ref['WLS_SIMULTANEOUS_FIBER_6'],
        err_msg='WLS_SIMULTANEOUS_FIBER_6',
    )
    np.testing.assert_array_equal(
        ad_out[0].REDUCED_ORDERS_FIBER_6, ref['REDUCED_ORDERS_FIBER_6'],
        err_msg='REDUCED_ORDERS_FIBER_6',
    )


# -- Recipes to create the inputs ----------------------------------------------
def _module_data_dir():
    root = os.environ.get('DRAGONS_TEST')
    if root is None:
        raise RuntimeError('DRAGONS_TEST environment variable not set')
    module_name = os.path.splitext(os.path.basename(__file__))[0]
    return os.path.join(
        root, 'maroonxdr', 'maroonx', 'echelle_extraction', module_name
    )


def create_inputs_recipe():
    """Copy the blessed reduced science products into inputs/."""
    import shutil
    from pathlib import Path

    inputs_dir = Path(_module_data_dir()) / 'inputs'
    inputs_dir.mkdir(parents=True, exist_ok=True)
    preprocessed_dir = Path(os.environ['DRAGONS_TEST']) / 'preprocessed_files'

    for filename in datasets:
        shutil.copy2(preprocessed_dir / filename, inputs_dir / filename)
        print(f'Copied to inputs/:\n    {filename}')


if __name__ == '__main__':
    import sys

    if '--create-inputs' in sys.argv[1:]:
        create_inputs_recipe()
    else:
        pytest.main([__file__])
