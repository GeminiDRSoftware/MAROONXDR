"""Regression tests for the fitAndApplyEtalonWls primitive.

The input files are the wavecal products written by preprocess/wavecal.py.
Each one carries both the primitive's input state (PEAKS/POLY tables and
WLS_STATIC_* extensions) and its blessed output (the WLS_DYNAMIC_* arrays,
the PEAK_DATA table, and the DRIFT_FIBER_* header keywords). The test re-runs
the dynamic fit for all etalon fibers and compares against the stored output.

Usage:
    python -m maroonxdr.maroonx.tests.echelle_extraction.test_fit_and_apply_etalon_wls --create-inputs
"""

import os
from copy import deepcopy

import numpy as np
import pytest
from gempy.utils import logutils

import astrodata
import maroonx_instruments  # noqa - registers the MaroonX AstroData class
from maroonxdr.maroonx.primitives_maroonx_spectrum import MaroonXSpectrum

# Wavecal products, per arm.
datasets = [
    '20250717T163124Z_DEEEE_b_0010_wavecal.fits',
    '20250717T163124Z_DEEEE_r_0004_wavecal.fits',
]

# Etalon fibers of the DEEEE frames
FIBERS = [2, 3, 4, 5]

# PEAK_DATA columns to compare
PD_INT_COLS = ['FIBER', 'ORDER', 'M']
PD_FLOAT_COLS = ['OFFSET', 'CENTER', 'AMPLITUDE', 'SIGMA1', 'SIGMA2',
                        'WIDTH', 'WAVELENGTH_BY_THAR', 'DISPERSION_MPS',
                        'M_FRACTION', 'WAVELENGTH']


# -- Tests ---------------------------------------------------------------------
@pytest.mark.preprocessed_data
@pytest.mark.regression
@pytest.mark.parametrize('filename', datasets)
def test_fit_and_apply_etalon_wls(filename, path_to_inputs, change_working_dir):
    """Re-fitted dynamic wavelength solution is stable against reference."""
    ad = astrodata.open(os.path.join(path_to_inputs, filename))

    # Stored output from the reference run
    ref_wls = {fiber: getattr(ad[0], f'WLS_DYNAMIC_FIBER_{fiber}')
               for fiber in FIBERS}
    ref_peak_data = ad[0].PEAK_DATA
    ref_drifts = {fiber: ad[0].hdr[f'DRIFT_FIBER_{fiber}'] for fiber in FIBERS}

    with change_working_dir():
        logutils.config(file_name=f'log_{filename.replace(".fits", "")}.txt')
        ad_out = MaroonXSpectrum([]).fitAndApplyEtalonWls(
            [deepcopy(ad)], report=False
        )[0]

    # Dynamic wavelength arrays
    for fiber in FIBERS:
        wls = getattr(ad_out[0], f'WLS_DYNAMIC_FIBER_{fiber}')
        assert wls.shape == ref_wls[fiber].shape
        np.testing.assert_allclose(
            wls, ref_wls[fiber], rtol=1e-10, err_msg=f'WLS_DYNAMIC_FIBER_{fiber}'
        )

    # Updated peak table
    peak_data = ad_out[0].PEAK_DATA
    assert len(peak_data) == len(ref_peak_data)
    for col in PD_INT_COLS:
        np.testing.assert_array_equal(
            peak_data[col], ref_peak_data[col], err_msg=col
        )
    for col in PD_FLOAT_COLS:
        np.testing.assert_allclose(
            peak_data[col], ref_peak_data[col], rtol=1e-6, atol=1e-8, err_msg=col
        )

    # Measured drifts
    for fiber in FIBERS:
        assert ad_out[0].hdr[f'DRIFT_FIBER_{fiber}'] == ref_drifts[fiber]


# -- Recipe to create the inputs -----------------------------------------------
def create_inputs_recipe():
    """Copy the wavecal products into inputs/."""
    import shutil
    from pathlib import Path

    root = os.environ.get('DRAGONS_TEST')
    if root is None:
        raise RuntimeError('DRAGONS_TEST environment variable not set')

    preprocessed_dir = Path(root) / 'preprocessed_files'
    module_name = os.path.splitext(os.path.basename(__file__))[0]
    inputs_dir = (
        Path(root) / 'maroonxdr' / 'maroonx' / 'echelle_extraction'
        / module_name / 'inputs'
    )
    inputs_dir.mkdir(parents=True, exist_ok=True)

    for filename in datasets:
        shutil.copy2(preprocessed_dir / filename, inputs_dir / filename)
        print(f'Copied to inputs/:\n    {filename}')


if __name__ == '__main__':
    import sys

    if '--create-inputs' in sys.argv[1:]:
        create_inputs_recipe()
    else:
        pytest.main([__file__])
