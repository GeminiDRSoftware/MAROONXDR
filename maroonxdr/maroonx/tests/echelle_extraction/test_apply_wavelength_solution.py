"""Regression tests for the applyWavelengthSolution primitive.

The input files are the reduced science products written by
preprocess/science.py, together with its dynamic wavecal calibration. 
The reduced product carries both the primitive's input state
(PEAKS/POLY from the science etalon fiber, WLS_STATIC_* and the box-extracted
spectra) and its blessed output (the WLS_SIMULTANEOUS_* arrays and the
INSTRUME_DRIFT/RELATIVE_DRIFT keywords), which no later primitive in the
recipe modifies. The test re-runs the drift correction and compares against
the stored output.

Usage:
    python -m maroonxdr.maroonx.tests.echelle_extraction.test_apply_wavelength_solution --create-inputs
"""

import os
from copy import deepcopy

import numpy as np
import pytest
from gempy.utils import logutils

import astrodata
import maroonx_instruments  # noqa - registers the MaroonX AstroData class
from maroonxdr.maroonx.primitives_maroonx_spectrum import MaroonXSpectrum

# Reduced science product and its dynamic wavecal calibration, per arm.
datasets = [
    ('20250717T144308Z_SOOOE_b_0300_reduced.fits',
     '20250717T163124Z_DEEEE_b_0010_wavecal.fits'),
    ('20250717T144308Z_SOOOE_r_0300_reduced.fits',
     '20250717T163124Z_DEEEE_r_0004_wavecal.fits'),
]

# Science fibers and simultaneous-etalon reference fiber, as in the recipe.
FIBERS = [2, 3, 4]
REF_FIBER = 5


# -- Tests ---------------------------------------------------------------------
@pytest.mark.preprocessed_data
@pytest.mark.regression
@pytest.mark.parametrize('filename, wavecal', datasets)
def test_apply_wavelength_solution(
    filename, wavecal, path_to_inputs, change_working_dir
):
    """Run drift-corrected wavelength solution against reference."""
    ad = astrodata.open(os.path.join(path_to_inputs, filename))
    wavecal_path = os.path.join(path_to_inputs, wavecal)

    # Stored output from the reference run
    ref_wls = {fiber: getattr(ad[0], f'WLS_SIMULTANEOUS_FIBER_{fiber}')
               for fiber in FIBERS + [REF_FIBER]}
    ref_inst_drift = ad[0].hdr['INSTRUME_DRIFT']
    ref_rel_drift = ad[0].hdr['RELATIVE_DRIFT']

    with change_working_dir():
        logutils.config(file_name=f'log_{filename.replace(".fits", "")}.txt')
        ad_out = MaroonXSpectrum([]).applyWavelengthSolution(
            [deepcopy(ad)],
            wavecal=wavecal_path,
            fibers=FIBERS,
            ref_fiber=REF_FIBER,
            report=False,
        )[0]

    # Drift-corrected wavelength arrays, per fiber
    for fiber in FIBERS + [REF_FIBER]:
        wls = getattr(ad_out[0], f'WLS_SIMULTANEOUS_FIBER_{fiber}')
        assert wls.shape == ref_wls[fiber].shape
        np.testing.assert_allclose(
            wls, ref_wls[fiber], rtol=1e-10,
            err_msg=f'WLS_SIMULTANEOUS_FIBER_{fiber}',
        )

    # Measured drifts
    assert ad_out[0].hdr['INSTRUME_DRIFT'] == ref_inst_drift
    assert ad_out[0].hdr['RELATIVE_DRIFT'] == ref_rel_drift


# -- Recipe to create the inputs -----------------------------------------------
def create_inputs_recipe():
    """Copy the reduced science products and their wavecals into inputs/."""
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

    for names in datasets:
        for filename in names:
            shutil.copy2(preprocessed_dir / filename, inputs_dir / filename)
            print(f'Copied to inputs/:\n    {filename}')


if __name__ == '__main__':
    import sys

    if '--create-inputs' in sys.argv[1:]:
        create_inputs_recipe()
    else:
        pytest.main([__file__])
