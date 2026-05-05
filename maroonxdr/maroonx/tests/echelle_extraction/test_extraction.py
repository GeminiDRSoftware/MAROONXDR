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
        '20241124T041907Z_SOOOE_r_0300_reduced.fits',
        '20241124T041907Z_SOOOE_b_0300_reduced.fits',
    ],
)
def test_optimal_extracting_science_data(caplog, path_to_inputs, filename):
    """
    Test that new optimal extraction and box extraction is equal to
    a previously reduced extraction of the same intial dataset
    for a standard science frame  (aka 3 science fibers, sim cal and sky fiber).
    As is standard for the echelle extraction, all extensions should exist even
    if not called to be populated by data
    (e.g. optimal extraction of simcal fiber).

    Parameters
    ----------
    caplog : fixture
    filename : str
    """
    caplog.set_level(logging.DEBUG)

    # Resolve calibration paths from preprocessed_files/calibrations/
    arm = 'r' if '_r_' in filename else 'b'
    calib_root = Path(os.environ['DRAGONS_TEST']) / 'preprocessed_files' / 'calibrations'
    flat_path = str(calib_root / 'processed_flat' / _CALIB_FILES[arm]['flat'])
    dark_path = str(calib_root / 'processed_dark' / _CALIB_FILES[arm]['dark'])

    # test that optimal extraction is equal to a previously reduced extraction
    ad = astrodata.open(os.path.join(path_to_inputs, filename))
    p = MAROONXEchelle([deepcopy(ad)])
    p.extractStripes(
        flat=flat_path,
        dark=dark_path,
        dark_subtraction_skip_fibers=[5],
        straylight_removal_fibers=[5],
    )
    adout = p.optimalExtraction(dark=dark_path)
    adtest = adout[0]

    assert len(caplog.records) > 0
    assert any('extracted' in r.message for r in caplog.records)

    attrs = [
        'REDUCED_ORDERS_FIBER_',
        'OPTIMAL_REDUCED_FIBER_',
        'OPTIMAL_REDUCED_VAR_',
        'BOX_REDUCED_FIBER_',
        'BOX_REDUCED_VAR_',
        'BPM_FIBER_',
    ]
    for f in range(1, 6):
        for attr in attrs:
            actual = getattr(adtest, f'{attr}{f}', None)
            expected = getattr(ad, f'{attr}{f}', None)
            if actual is None and expected is None:
                continue
            np.testing.assert_array_equal(
                actual, expected,
                err_msg=f'{attr}{f} mismatch',
            )

# helper function to compare stripe dicts
def _compare_stripes(stripes1, stripes2, label, 
    assert_func=np.testing.assert_array_equal, **kwargs):
    """Compare two STRIPES dicts (fiber -> order -> sparse matrix)."""
    assert sorted(stripes1.keys()) == sorted(stripes2.keys()), \
        f'{label}: fiber keys differ'
    for fiber in sorted(stripes1.keys()):
        assert sorted(stripes1[fiber].keys()) == sorted(stripes2[fiber].keys()), \
            f'{label} {fiber}: order keys differ'
        for order in sorted(stripes1[fiber].keys()):
            a = stripes1[fiber][order].toarray()
            b = stripes2[fiber][order].toarray()
            assert_func(
                a, b,
                err_msg=f'{label} {fiber}/{order}',
                **kwargs,
            )


@pytest.mark.slow
@pytest.mark.preprocessed_data
@pytest.mark.parametrize(
    'filename',
    [
        '20241124T041907Z_SOOOE_r_0300_reduced.fits',
        '20241124T041907Z_SOOOE_b_0300_reduced.fits',
    ],
)
def test_extraction_determinism_no_straylight(caplog, path_to_inputs, filename):
    """
    Test that extractStripes is deterministic without straylight removal.

    Runs extractStripes twice with no dark subtraction and no straylight
    removal, then compares the STRIPES output for bitwise equality.

    Parameters
    ----------
    caplog : fixture
    filename : str
    """
    caplog.set_level(logging.DEBUG)

    arm = 'r' if '_r_' in filename else 'b'
    calib_root = Path(os.environ['DRAGONS_TEST']) / 'preprocessed_files' / 'calibrations'
    flat_path = str(calib_root / 'processed_flat' / _CALIB_FILES[arm]['flat'])

    ad = astrodata.open(os.path.join(path_to_inputs, filename))

    # First run — no dark, no straylight
    p1 = MAROONXEchelle([deepcopy(ad)])
    p1.extractStripes(
        flat=flat_path,
        dark_subtraction_skip_fibers=[1, 2, 3, 4, 5],
        straylight_removal_fibers=[],
    )
    r1 = p1.streams['main'][0][0]

    # Second run — identical parameters
    p2 = MAROONXEchelle([deepcopy(ad)])
    p2.extractStripes(
        flat=flat_path,
        dark_subtraction_skip_fibers=[1, 2, 3, 4, 5],
        straylight_removal_fibers=[],
    )
    r2 = p2.streams['main'][0][0]

    _compare_stripes(r1.STRIPES, r2.STRIPES, 'STRIPES')
    _compare_stripes(r1.F_STRIPES, r2.F_STRIPES, 'F_STRIPES')
    _compare_stripes(r1.STRIPES_MASKS, r2.STRIPES_MASKS, 'STRIPES_MASKS')


@pytest.mark.slow
@pytest.mark.preprocessed_data
@pytest.mark.parametrize(
    'filename',
    [
        '20241124T041907Z_SOOOE_r_0300_reduced.fits',
        '20241124T041907Z_SOOOE_b_0300_reduced.fits',
    ],
)
def test_extraction_determinism_with_straylight(caplog, path_to_inputs, filename):
    """
    Test that extractStripes is deterministic with straylight removal enabled.

    Runs extractStripes twice with dark subtraction and straylight removal
    (same params as the science recipe), then compares the STRIPES output.
    If this test fails but test_extraction_determinism_no_straylight passes,
    that confirms Background2D straylight fitting is the source of
    non-determinism.

    Parameters
    ----------
    caplog : fixture
    filename : str
    """
    caplog.set_level(logging.DEBUG)

    arm = 'r' if '_r_' in filename else 'b'
    calib_root = Path(os.environ['DRAGONS_TEST']) / 'preprocessed_files' / 'calibrations'
    flat_path = str(calib_root / 'processed_flat' / _CALIB_FILES[arm]['flat'])
    dark_path = str(calib_root / 'processed_dark' / _CALIB_FILES[arm]['dark'])

    ad = astrodata.open(os.path.join(path_to_inputs, filename))

    # First run — with dark subtraction and straylight (science recipe params)
    p1 = MAROONXEchelle([deepcopy(ad)])
    p1.extractStripes(
        flat=flat_path,
        dark=dark_path,
        dark_subtraction_skip_fibers=[5],
        straylight_removal_fibers=[5],
    )
    r1 = p1.streams['main'][0][0]

    # Second run — identical parameters
    p2 = MAROONXEchelle([deepcopy(ad)])
    p2.extractStripes(
        flat=flat_path,
        dark=dark_path,
        dark_subtraction_skip_fibers=[5],
        straylight_removal_fibers=[5],
    )
    r2 = p2.streams['main'][0][0]

    _compare_stripes(r1.STRIPES, r2.STRIPES, 'STRIPES')
    _compare_stripes(r1.F_STRIPES, r2.F_STRIPES, 'F_STRIPES')
    _compare_stripes(r1.STRIPES_MASKS, r2.STRIPES_MASKS, 'STRIPES_MASKS')


# Faster tests
@pytest.mark.preprocessed_data
def test_boxExtraction(ad_echelle):
    """Box-extracted flux equals column sum of the stripe."""
    signal = 100.0
    slit_height = 10
    expected_flux = signal * 2 * slit_height

    p = MAROONXEchelle([])
    result = p.boxExtraction([ad_echelle])[0]

    for fiber in range(1, 6):
        box = getattr(result[0], f'BOX_REDUCED_FIBER_{fiber}')
        assert box.shape[0] > 0
        np.testing.assert_allclose(box, expected_flux, atol=1e-4,
                                   err_msg=f"fiber {fiber}")


@pytest.mark.preprocessed_data
def test_optimalExtraction(ad_echelle):
    """Optimal extraction recovers flux from uniform stripes."""
    signal = 100.0
    slit_height = 10
    expected_box_flux = signal * 2 * slit_height

    p = MAROONXEchelle([])
    result = p.optimalExtraction([ad_echelle], 
        back_var=0.0,
        optimal_extraction_fibers=[3])[0]

    # Fiber 3: optimally extracted
    opt = getattr(result[0], 'OPTIMAL_REDUCED_FIBER_3')
    box = getattr(result[0], 'BOX_REDUCED_FIBER_3')
    assert opt.shape != (1, 1)
    assert box.shape != (1, 1)
    np.testing.assert_allclose(box, expected_box_flux, atol=1e-4,
                               err_msg="box fiber 3")
    np.testing.assert_allclose(opt, expected_box_flux, rtol=1e-4,
                               err_msg="optimal fiber 3")

    # Other fibers: no optimal extraction, only placeholder
    # box is always extracted
    for fiber in [1, 2, 4, 5]:
        opt = getattr(result[0], f'OPTIMAL_REDUCED_FIBER_{fiber}')
        box = getattr(result[0], f'BOX_REDUCED_FIBER_{fiber}')
        assert opt.shape == (1, 1), f"fiber {fiber} should be placeholder"
        assert box.shape != (1, 1), f"fiber {fiber} should have box extraction"


if __name__ == '__main__':
    pytest.main()
