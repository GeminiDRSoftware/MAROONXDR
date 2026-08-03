"""Regression tests for the getPeaksAndPolynomials primitive.

The input files are the wavecal products written by preprocess/wavecal.py.
Each one carries both the primitive's input state (the BOX_REDUCED_*
extensions from boxExtraction) and its blessed output (the PEAKS and POLY
tables, which no later primitive in the recipe modifies). The test re-fits a
small fiber/order subset in serial mode and compares against the stored
tables; a full re-fit takes minutes per arm and is not needed since the fit
is independent per fiber/order.

Usage:
    python -m maroonxdr.maroonx.tests.echelle_extraction.test_get_peaks_and_polynomials --create-inputs
"""

import os
from copy import deepcopy

import numpy as np
import pytest
from gempy.utils import logutils

import astrodata
import maroonx_instruments  # noqa - registers the MaroonX AstroData class
from maroonxdr.maroonx.primitives_maroonx_spectrum import MaroonXSpectrum

# Wavecal products and select the orders to re-fit
datasets = [
    ('20250717T163124Z_DEEEE_b_0010_wavecal.fits', [100, 108]),
    ('20250717T163124Z_DEEEE_r_0004_wavecal.fits', [75, 81]),
]

FIBER = 5

# PEAKS/POLY columns compared exactly (identifiers, fit status) and to
# within floating-point tolerance (fitted quantities).
PEAKS_INT_COLS = ['FIBER', 'ORDER']
PEAKS_FLOAT_COLS = ['OFFSET', 'CENTER', 'AMPLITUDE', 'SIGMA1', 'SIGMA2',
                    'WIDTH']
POLY_INT_COLS = ['fiber', 'order']
POLY_FLOAT_COLS = ['offset', 'fitrange', 'sigma_l_coefficients',
                   'sigma_r_coefficients', 'width_coefficients']


# -- Tests ---------------------------------------------------------------------
@pytest.mark.preprocessed_data
@pytest.mark.regression
@pytest.mark.parametrize('filename, orders', datasets)
def test_get_peaks_and_polynomials(
    filename, orders, path_to_inputs, change_working_dir
):
    """Re-fitted etalon peaks are stable against the blessed tables."""
    ad = astrodata.open(os.path.join(path_to_inputs, filename))

    # Stored tables from the blessed run, restricted to the re-fitted subset
    ref_peaks = ad[0].PEAKS
    ref_peaks = ref_peaks[
        (ref_peaks['FIBER'] == FIBER) & np.isin(ref_peaks['ORDER'], orders)
    ]
    ref_poly = ad[0].POLY
    ref_poly = ref_poly[
        (ref_poly['fiber'] == FIBER) & np.isin(ref_poly['order'], orders)
    ]

    with change_working_dir():
        logutils.config(file_name=f'log_{filename.replace(".fits", "")}.txt')
        ad_out = MaroonXSpectrum([]).getPeaksAndPolynomials(
            [deepcopy(ad)], fibers=[FIBER], orders=orders, multithreading=False
        )[0]

    peaks = ad_out[0].PEAKS
    poly = ad_out[0].POLY

    assert len(peaks) == len(ref_peaks)
    for col in PEAKS_INT_COLS:
        np.testing.assert_array_equal(peaks[col], ref_peaks[col], err_msg=col)
    for col in PEAKS_FLOAT_COLS:
        np.testing.assert_allclose(
            peaks[col], ref_peaks[col], rtol=1e-6, atol=1e-8, err_msg=col
        )

    assert len(poly) == len(ref_poly)
    for col in POLY_INT_COLS:
        np.testing.assert_array_equal(poly[col], ref_poly[col], err_msg=col)
    for col in POLY_FLOAT_COLS:
        np.testing.assert_allclose(
            poly[col], ref_poly[col], rtol=1e-6, atol=1e-8, err_msg=col
        )


# -- Recipe to create the inputs -----------------------------------------------
def create_inputs_recipe():
    """Copy the wavecal products into inputs/.

    The products written by preprocess/wavecal.py already contain both the
    input state and the blessed PEAKS/POLY tables, so staging is a copy from
    ``$DRAGONS_TEST/preprocessed_files/``.
    """
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

    for filename, _ in datasets:
        shutil.copy2(preprocessed_dir / filename, inputs_dir / filename)
        print(f'Copied to inputs/:\n    {filename}')


if __name__ == '__main__':
    import sys

    if '--create-inputs' in sys.argv[1:]:
        create_inputs_recipe()
    else:
        pytest.main([__file__])
