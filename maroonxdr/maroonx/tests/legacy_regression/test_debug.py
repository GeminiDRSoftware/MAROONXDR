"""
Debug tests to investigate precision differences between
new (DRAGONS) and legacy box extraction pipelines.

Hypothesis: the flat box extraction differs because the new pipeline
applies the science frame's BPM when box-summing the flat stripes,
while the legacy pipeline used the flat's own mask.

These tests isolate each component of the load_recordings division:
  data = box_extraction / flat_box_extraction

to identify exactly where the discrepancy originates.
"""

import astrodata
import h5py
import numpy as np
import pytest
import tables
from gempy.utils import logutils

import maroonx_instruments  # noqa : important to load adclass tags
from maroonxdr.maroonx.primitives_maroonx_spectrum import MaroonXSpectrum

from maroonxdr.maroonx.primitives_maroonx_echelle import (
    _extract_single_stripe,
    _box_extract_single_stripe,
)
from . import legacy_adapter

_CALIB_FILES = {
    'BLUE': {
        'flat': 'processed_flat/20241114T190714Z_DDDDF_b_0007_DFFFF_flat.fits',
        'dark': 'processed_dark/20241124T041907Z_SOOOE_b_0300_synth_dark.fits',
    },
}

logutils.config(file_name="test_debug.log", mode="debug", stomp=True)
log = logutils.get_logger("test_debug")
log.setLevel("DEBUG")

FIBER = 2
ORDER = '111'
ETALON = '20241124T030227Z_DEEEE_b_0030'


@pytest.mark.preprocessed_data
def test_flat_box_extraction_comparison(
    path_to_legacy_flats, preprocessed_files_path
):
    """
    Test 1: Compare the flat box extraction values directly.

    If the flat box extractions differ, then the division
    (science / flat) will differ even if the science numerators are identical.
    """
    arm = 'BLUE'
    calib_root = preprocessed_files_path / 'calibrations'
    flat_path = str(calib_root / _CALIB_FILES[arm]['flat'])
    dark_path = str(calib_root / _CALIB_FILES[arm]['dark'])

    raw_file = preprocessed_files_path / (ETALON + '.fits')
    adinput = [astrodata.open(raw_file)]

    p = MaroonXSpectrum(adinput)
    p.prepare()
    p.checkArm()
    p.addDQ()
    p.overscanCorrect()
    p.correctImageOrientation()
    p.addVAR(read_noise=True, poisson_noise=True)

    p.extractStripes(flat=flat_path, dark=dark_path)
    adout = p.boxExtraction()
    ad = adout[0]

    # New flat box extraction for fiber/order
    new_orders = ad[0].REDUCED_ORDERS_FIBER_2.astype(int)
    order_idx = list(new_orders).index(int(ORDER))
    new_flat_box = ad[0].BOX_REDUCED_FLAT_2[order_idx]

    # Legacy flat box extraction (from the flat's own HDF)
    old_flat_file = (
        path_to_legacy_flats
        / '20241114T19_masterflat_backgroundsubtracted_FFFFF_b_0007.hdf'
    )
    with tables.open_file(str(old_flat_file), 'r') as h5flat:
        legacy_flat_box = np.array(
            h5flat.get_node(f'/box_extraction/fiber_{FIBER}/{ORDER}')
        )

    # Handle NaN values: legacy may have NaN at detector edges
    new_nan = np.isnan(new_flat_box)
    leg_nan = np.isnan(legacy_flat_box)
    both_valid = ~new_nan & ~leg_nan

    print(f"\n--- Flat box extraction comparison (fiber {FIBER}, order {ORDER}) ---")
    print(f"  Shape: new={new_flat_box.shape}, legacy={legacy_flat_box.shape}")
    print(f"  NaN pixels: new={new_nan.sum()}, legacy={leg_nan.sum()}")
    print(f"  NaN only in new: {(new_nan & ~leg_nan).sum()}")
    print(f"  NaN only in legacy: {(~new_nan & leg_nan).sum()}")

    diff = np.abs(new_flat_box[both_valid] - legacy_flat_box[both_valid])
    n_diff = np.sum(diff > 1e-10)
    max_diff = np.max(diff) if len(diff) > 0 else 0.0
    print(f"  Valid pixels compared: {both_valid.sum()} / {len(new_flat_box)}")
    print(f"  Pixels differing (>1e-10): {n_diff}")
    print(f"  Max absolute difference: {max_diff:.2e}")
    print(f"  New flat range: [{np.nanmin(new_flat_box):.4f}, {np.nanmax(new_flat_box):.4f}]")
    print(f"  Legacy flat range: [{np.nanmin(legacy_flat_box):.4f}, {np.nanmax(legacy_flat_box):.4f}]")

    # This assertion is expected to FAIL if the hypothesis is correct
    np.testing.assert_allclose(
        new_flat_box[both_valid], legacy_flat_box[both_valid], rtol=0, atol=1e-10,
        err_msg="Flat box extractions differ — mask hypothesis likely correct"
    )


@pytest.mark.preprocessed_data
def test_science_box_extraction_comparison(
    path_to_legacy_wavecal, preprocessed_files_path
):
    """
    Test 2: Compare the science (wavecal) box extraction values directly.

    This isolates the numerator of the division. If only the flat differs
    but the science box extraction matches, it confirms the mask hypothesis.
    """
    arm = 'BLUE'
    calib_root = preprocessed_files_path / 'calibrations'
    flat_path = str(calib_root / _CALIB_FILES[arm]['flat'])
    dark_path = str(calib_root / _CALIB_FILES[arm]['dark'])

    raw_file = preprocessed_files_path / (ETALON + '.fits')
    adinput = [astrodata.open(raw_file)]

    p = MaroonXSpectrum(adinput)
    p.prepare()
    p.checkArm()
    p.addDQ()
    p.overscanCorrect()
    p.correctImageOrientation()
    p.addVAR(read_noise=True, poisson_noise=True)

    p.extractStripes(flat=flat_path, dark=dark_path)
    adout = p.boxExtraction()
    ad = adout[0]

    # New science box extraction
    new_orders = ad[0].REDUCED_ORDERS_FIBER_2.astype(int)
    order_idx = list(new_orders).index(int(ORDER))
    new_science_box = ad[0].BOX_REDUCED_FIBER_2[order_idx]

    # Legacy science box extraction (from the wavecal HDF)
    old_file = path_to_legacy_wavecal / (ETALON + '.hdf')
    with tables.open_file(str(old_file), 'r') as h5f:
        legacy_science_box = np.array(
            h5f.get_node(f'/box_extraction/fiber_{FIBER}/{ORDER}')
        )

    # Both pipelines now produce NaN for fully-masked columns
    new_nan = np.isnan(new_science_box)
    leg_nan = np.isnan(legacy_science_box)
    both_valid = ~new_nan & ~leg_nan

    diff = np.abs(new_science_box[both_valid] - legacy_science_box[both_valid])
    n_diff = np.sum(diff > 1e-10)
    max_diff = np.max(diff) if len(diff) > 0 else 0.0
    print(f"\n--- Science box extraction comparison (fiber {FIBER}, order {ORDER}) ---")
    print(f"  Shape: new={new_science_box.shape}, legacy={legacy_science_box.shape}")
    print(f"  NaN pixels: new={new_nan.sum()}, legacy={leg_nan.sum()}")
    print(f"  NaN agreement: {(new_nan == leg_nan).all()}")
    print(f"  Valid pixels compared: {both_valid.sum()} / {len(new_science_box)}")
    print(f"  Pixels differing (>1e-10): {n_diff} / {both_valid.sum()}")
    print(f"  Max absolute difference: {max_diff:.2e}")

    # DRAGONS may have extra NaN pixels from the zero→NaN fix in boxExtraction.
    # Legacy never had that conversion, so new_nan can be a superset of leg_nan.
    # What matters: valid pixels must match exactly.
    assert (~leg_nan).all() or (new_nan[leg_nan]).all(), (
        "Legacy has NaN where DRAGONS does not — unexpected"
    )
    np.testing.assert_allclose(
        new_science_box[both_valid], legacy_science_box[both_valid],
        rtol=0, atol=1e-10,
        err_msg="Science box extractions differ on valid pixels"
    )


@pytest.mark.preprocessed_data
def test_bpm_comparison(
    legacy_test_root, preprocessed_files_path
):
    """
    Test 3: Compare the BPM used in box extraction.

    The new pipeline applies the science frame's DQ/mask to both science
    and flat stripes during box extraction. The legacy pipeline box-extracted
    the flat independently with a static config BPM (from config_b.hdf).

    This test extracts the stripe-level mask from both sources and compares
    them pixel-by-pixel to identify exactly where they differ.
    """
    arm = 'BLUE'
    calib_root = preprocessed_files_path / 'calibrations'
    flat_path = str(calib_root / _CALIB_FILES[arm]['flat'])
    dark_path = str(calib_root / _CALIB_FILES[arm]['dark'])

    # --- Load legacy config BPM ---
    config_path = (
        legacy_test_root
        / 'MaroonX_spectra_reduced'
        / 'Maroonx_configfiles'
        / '202411xx'
        / 'config_b.hdf'
    )
    with h5py.File(str(config_path), 'r') as cf:
        legacy_bpm_raw = np.array(cf['bad_pixel_map'])   # (4400, 4400), 1=good
        valid = np.array(cf['valid'])                     # overscan map

    # Remove overscan (same as legacy utils.Configuration.remove_overscan)
    legacy_bpm = legacy_bpm_raw[valid == 1].reshape(
        (np.sum(valid[:, 2000]), np.sum(valid[2000]))
    )
    # Apply orientation correction for blue arm (flip both axes)
    legacy_bpm = np.flip(legacy_bpm, (0, 1))

    # --- Run new pipeline up to extractStripes ---
    raw_file = preprocessed_files_path / (ETALON + '.fits')
    adinput = [astrodata.open(raw_file)]

    p = MaroonXSpectrum(adinput)
    p.prepare()
    p.checkArm()
    p.addDQ()
    p.overscanCorrect()
    p.correctImageOrientation()
    p.addVAR(read_noise=True, poisson_noise=True)

    p.extractStripes(flat=flat_path, dark=dark_path)

    ad = p.streams['main'][0]
    stripes_masks = ad[0].STRIPES_MASKS

    # Reconstruct p_id dict from flat's STRIPES_ID table
    # (same logic as extractStripes lines 434-442)
    flat_ad = astrodata.open(flat_path)
    p_id_table = flat_ad[0].STRIPES_ID
    p_id = {}
    for i, f in enumerate(flat_ad[0].STRIPES_FIBERS):
        p_id[f'fiber_{f}'] = dict(
            (colname, p_id_table[6 * i : 6 * (i + 1)][colname].data)
            for colname in p_id_table.colnames
        )

    slit_height = 10  # default from extractStripes parameters
    fiber_key = f'fiber_{FIBER}'

    print(f"\n--- BPM comparison (fiber {FIBER}, order {ORDER}) ---")
    print(f"  Legacy config BPM shape (after overscan+orientation): {legacy_bpm.shape}")
    print(f"  Legacy config BPM bad pixels: {np.sum(legacy_bpm == 0)}")

    # Get the polynomial for this fiber/order
    poly_coeffs = p_id[fiber_key][ORDER]

    # Extract stripe-level mask from legacy config BPM
    legacy_stripe_mask = _extract_single_stripe(
        legacy_bpm.astype(int), poly_coeffs, slit_height
    )

    # Science DQ stripe mask (from the new pipeline)
    science_stripe_mask = stripes_masks[fiber_key][ORDER]

    # Convert sparse to dense for comparison
    legacy_dense = np.array(legacy_stripe_mask.toarray())
    science_dense = np.array(science_stripe_mask.toarray())

    print(f"  Stripe mask shape: {legacy_dense.shape}")
    print(f"  Legacy stripe good pixels: {legacy_dense.sum()}")
    print(f"  Science stripe good pixels: {science_dense.sum()}")

    # Find pixels that differ
    diff_mask = legacy_dense != science_dense
    n_diff_pixels = np.sum(diff_mask)
    print(f"  Pixels where masks differ: {n_diff_pixels}")

    if n_diff_pixels > 0:
        diff_rows, diff_cols = np.where(diff_mask)
        # Columns that differ (dispersion direction)
        diff_columns = np.unique(diff_cols)
        print(f"  Dispersion columns affected: {len(diff_columns)}")
        print(f"  Affected column indices: {diff_columns[:20]}{'...' if len(diff_columns) > 20 else ''}")

        # For each affected column, show the mask values
        for col in diff_columns[:5]:
            leg_col = legacy_dense[:, col]
            sci_col = science_dense[:, col]
            print(f"    Column {col}: legacy_good={leg_col.sum()}, "
                  f"science_good={sci_col.sum()}, "
                  f"diff={leg_col.sum() - sci_col.sum()}")

        # Show impact on box extraction: the column sum difference
        legacy_colsum = np.array(legacy_stripe_mask.sum(axis=0)).flatten()
        science_colsum = np.array(science_stripe_mask.sum(axis=0)).flatten()
        colsum_diff = legacy_colsum - science_colsum
        print(f"\n  Column sum differences (legacy - science):")
        print(f"    Non-zero columns: {np.sum(colsum_diff != 0)}")
        print(f"    Max diff: {colsum_diff.max()}, Min diff: {colsum_diff.min()}")

        print("\n  -> Masks DIFFER — this explains flat box extraction discrepancy")
    else:
        print("\n  -> Masks are IDENTICAL for this fiber/order")
        print("     Mask hypothesis is NOT the cause for this order")

    # Also compare across all orders for this fiber
    print(f"\n--- Summary across all orders for fiber {FIBER} ---")
    total_diff_cols = 0
    for order_key, poly in p_id[fiber_key].items():
        leg_stripe = _extract_single_stripe(
            legacy_bpm.astype(int), poly, slit_height
        )
        sci_stripe = stripes_masks[fiber_key][order_key]
        leg_arr = np.array(leg_stripe.toarray())
        sci_arr = np.array(sci_stripe.toarray())
        diff = leg_arr != sci_arr
        n_diff = np.sum(diff)
        if n_diff > 0:
            n_cols = len(np.unique(np.where(diff)[1]))
            total_diff_cols += n_cols
            print(f"  Order {order_key}: {n_diff} pixels differ across {n_cols} columns")
    if total_diff_cols == 0:
        print("  All orders: masks identical")
    else:
        print(f"  Total columns with mask differences: {total_diff_cols}")


@pytest.mark.preprocessed_data
def test_load_recordings_decomposed(
    path_to_legacy_wavecal, path_to_legacy_flats, preprocessed_files_path
):
    """
    Test 4: Decompose the load_recordings comparison into numerator and denominator.

    For the specific fiber/order that fails, show:
    - Whether science box extractions match
    - Whether flat box extractions match
    - How much the final ratio differs

    This pinpoints whether the problem is in the numerator, denominator, or both.
    """
    arm = 'BLUE'
    calib_root = preprocessed_files_path / 'calibrations'
    flat_path = str(calib_root / _CALIB_FILES[arm]['flat'])
    dark_path = str(calib_root / _CALIB_FILES[arm]['dark'])

    raw_file = preprocessed_files_path / (ETALON + '.fits')
    adinput = [astrodata.open(raw_file)]

    p = MaroonXSpectrum(adinput)
    p.prepare()
    p.checkArm()
    p.addDQ()
    p.overscanCorrect()
    p.correctImageOrientation()
    p.addVAR(read_noise=True, poisson_noise=True)

    p.extractStripes(flat=flat_path, dark=dark_path)
    adout = p.boxExtraction()
    ad = adout[0]

    new_orders = ad[0].REDUCED_ORDERS_FIBER_2.astype(int)

    # Legacy data
    old_file = path_to_legacy_wavecal / (ETALON + '.hdf')
    old_flat_file = (
        path_to_legacy_flats
        / '20241114T19_masterflat_backgroundsubtracted_FFFFF_b_0007.hdf'
    )

    print(f"\n--- Decomposed comparison for all orders of fiber {FIBER} ---")
    print(f"{'Order':>6} | {'Sci diff':>10} | {'Flat diff':>10} | {'Ratio diff':>10} | {'Ratio max':>10}")
    print("-" * 65)

    for order_idx, order in enumerate(new_orders):
        order_str = str(order)

        new_sci = ad[0].BOX_REDUCED_FIBER_2[order_idx]
        new_flat = ad[0].BOX_REDUCED_FLAT_2[order_idx]
        new_ratio = new_sci / new_flat

        with tables.open_file(str(old_file), 'r') as h5f:
            legacy_sci = np.array(
                h5f.get_node(f'/box_extraction/fiber_{FIBER}/{order_str}')
            )
        with tables.open_file(str(old_flat_file), 'r') as h5flat:
            legacy_flat = np.array(
                h5flat.get_node(f'/box_extraction/fiber_{FIBER}/{order_str}')
            )
        legacy_ratio = legacy_sci / legacy_flat

        sci_maxdiff = np.max(np.abs(new_sci - legacy_sci))
        flat_maxdiff = np.max(np.abs(new_flat - legacy_flat))
        ratio_maxdiff = np.max(np.abs(new_ratio - legacy_ratio))
        ratio_n_fail = np.sum(~np.isclose(new_ratio, legacy_ratio, rtol=0, atol=1e-5))

        print(
            f"{order_str:>6} | {sci_maxdiff:>10.2e} | {flat_maxdiff:>10.2e} | "
            f"{ratio_maxdiff:>10.2e} | fails={ratio_n_fail}"
        )


# =====================================================================
# Tests 5-8: Investigate the source of the flat box extraction
# discrepancy between DRAGONS and legacy pipelines.
# =====================================================================

# Paths to the legacy FFFFF flat files (relative to path_to_legacy_flats)
LEGACY_FFFFF_FITS = (
    '20241114T19_masterflat_backgroundsubtracted_FFFFF_b_0007.fits'
)
LEGACY_FFFFF_HDF = (
    '20241114T19_masterflat_backgroundsubtracted_FFFFF_b_0007.hdf'
)


@pytest.mark.preprocessed_data
def test_2d_flat_image_comparison(
    path_to_legacy_flats, preprocessed_files_path
):
    """
    Test 5: Compare the 2D flat images directly (before any stripe extraction).

    This tests whether the difference is already present in the 2D data.
    Loads the DRAGONS processed flat FITS and the legacy FFFFF background-
    subtracted FITS and compares pixel values.

    Finding: The two images differ at float32 precision (~90% of pixels,
    max abs diff ~2.8). The discrepancy originates in the flat processing
    pipeline (overscan, straylight removal, or flat combination) — not in
    stripe extraction or box summation.
    """
    from astropy.io import fits

    arm = 'BLUE'
    calib_root = preprocessed_files_path / 'calibrations'
    dragons_flat_path = str(calib_root / _CALIB_FILES[arm]['flat'])
    legacy_flat_path = str(path_to_legacy_flats / LEGACY_FFFFF_FITS)

    # Load DRAGONS flat (float64)
    dragons_ad = astrodata.open(dragons_flat_path)
    dragons_data = dragons_ad[0].data

    # Load legacy FFFFF flat (float32)
    with fits.open(legacy_flat_path) as hdul:
        legacy_data_f32 = hdul[0].data.copy()
    legacy_data_f64 = legacy_data_f32.astype(np.float64)

    print(f"\n--- Test 5: 2D flat image comparison ---")
    print(f"  DRAGONS flat dtype: {dragons_data.dtype}, shape: {dragons_data.shape}")
    print(f"  Legacy flat dtype:  {legacy_data_f32.dtype}, shape: {legacy_data_f32.shape}")

    # Compare at float64 precision (shows the dtype-induced difference)
    diff_f64 = np.abs(dragons_data - legacy_data_f64)
    n_diff_f64 = np.sum(diff_f64 > 1e-10)
    print(f"\n  Float64 comparison:")
    print(f"    Max abs diff: {diff_f64.max():.6e}")
    print(f"    Pixels differing (>1e-10): {n_diff_f64} / {diff_f64.size}")

    if n_diff_f64 > 0:
        rdiff = diff_f64[diff_f64 > 1e-10] / np.abs(legacy_data_f64[diff_f64 > 1e-10])
        print(f"    Max relative diff: {rdiff.max():.6e}")

    # Compare at float32 precision (the ground truth comparison)
    dragons_data_f32 = dragons_data.astype(np.float32)
    diff_f32 = np.abs(dragons_data_f32 - legacy_data_f32)
    n_diff_f32 = np.sum(diff_f32 > 0)

    print(f"\n  Float32 comparison (both cast to float32):")
    print(f"    Max abs diff: {diff_f32.max():.6e}")
    print(f"    Pixels differing (>0): {n_diff_f32} / {diff_f32.size}")

    if n_diff_f32 == 0:
        print("\n  -> CONCLUSION: 2D flat images are IDENTICAL at float32 precision.")
        print("     The flat box extraction discrepancy is caused entirely by")
        print("     float32 vs float64 dtype difference, not by any processing")
        print("     difference (overscan, straylight, flat combination, etc.).")

    # The images must be identical at float32 precision
    np.testing.assert_array_equal(
        dragons_data_f32, legacy_data_f32,
        err_msg="2D flat images differ even at float32 precision"
    )


@pytest.mark.preprocessed_data
def test_legacy_flat_reextraction(
    path_to_legacy_flats, preprocessed_files_path
):
    """
    Test 6: Box-extract the legacy FFFFF flat 2D image using new pipeline code.

    Takes the legacy FFFFF background-subtracted FITS, extracts a stripe using
    _extract_single_stripe, then box-sums with _box_extract_single_stripe.
    Compares the result to the legacy HDF5 /box_extraction/.

    If they match, the extraction code is correct and the only difference
    is in the input 2D data (confirmed by Test 5 as a dtype issue).
    """
    from astropy.io import fits

    legacy_flat_path = str(path_to_legacy_flats / LEGACY_FFFFF_FITS)
    legacy_hdf_path = str(path_to_legacy_flats / LEGACY_FFFFF_HDF)

    # Load legacy FFFFF 2D image (float32, as legacy used it)
    # Use native byte order for scipy.sparse compatibility
    with fits.open(legacy_flat_path) as hdul:
        legacy_2d = hdul[0].data.astype(np.float32)  # native byte order

    # Load legacy config BPM for masking
    config_path = (
        path_to_legacy_flats.parent.parent.parent
        / 'Maroonx_configfiles'
        / '202411xx'
        / 'config_b.hdf'
    )

    import h5py as h5
    with h5.File(str(config_path), 'r') as cf:
        legacy_bpm_raw = np.array(cf['bad_pixel_map'])
        valid = np.array(cf['valid'])

    # Remove overscan from BPM (same as legacy Configuration.remove_overscan)
    legacy_bpm = legacy_bpm_raw[valid == 1].reshape(
        (np.sum(valid[:, 2000]), np.sum(valid[2000]))
    )
    # Blue arm orientation correction
    legacy_bpm = np.flip(legacy_bpm, (0, 1))

    # Get polynomial coefficients from the DRAGONS flat
    calib_root = preprocessed_files_path / 'calibrations'
    flat_ad = astrodata.open(str(calib_root / _CALIB_FILES['BLUE']['flat']))
    p_id_table = flat_ad[0].STRIPES_ID
    p_id = {}
    for i, f in enumerate(flat_ad[0].STRIPES_FIBERS):
        p_id[f'fiber_{f}'] = dict(
            (colname, p_id_table[6 * i : 6 * (i + 1)][colname].data)
            for colname in p_id_table.colnames
        )

    slit_height = 10
    fiber_key = f'fiber_{FIBER}'
    poly_coeffs = p_id[fiber_key][ORDER]

    # Extract stripe from legacy 2D image using new code
    f_stripe = _extract_single_stripe(legacy_2d, poly_coeffs, slit_height)
    bpm_stripe = _extract_single_stripe(
        legacy_bpm.astype(int), poly_coeffs, slit_height
    )

    # Box extract using new code
    reextracted = _box_extract_single_stripe(f_stripe, bpm_stripe)

    # Load legacy box extraction from HDF5
    with tables.open_file(legacy_hdf_path, 'r') as h5flat:
        legacy_box = np.array(
            h5flat.get_node(f'/box_extraction/fiber_{FIBER}/{ORDER}')
        )

    # Handle NaN
    both_valid = ~np.isnan(reextracted) & ~np.isnan(legacy_box)
    diff = np.abs(reextracted[both_valid] - legacy_box[both_valid])

    print(f"\n--- Test 6: Re-extract legacy flat with new code ---")
    print(f"  Legacy 2D dtype: {legacy_2d.dtype}")
    print(f"  Re-extracted shape: {reextracted.shape}")
    print(f"  Legacy HDF box shape: {legacy_box.shape}")
    print(f"  Valid pixels: {both_valid.sum()}")
    print(f"  Max abs diff: {diff.max():.6e}")
    print(f"  Pixels differing (>1e-10): {np.sum(diff > 1e-10)}")

    if np.sum(diff > 1e-10) == 0:
        print("\n  -> CONCLUSION: New extraction code reproduces legacy box")
        print("     extraction exactly when given the same float32 input.")
    else:
        n_diff = np.sum(diff > 1e-10)
        print(f"\n  -> {n_diff} pixels differ. Extraction code may differ from legacy.")

    np.testing.assert_allclose(
        reextracted[both_valid], legacy_box[both_valid], rtol=0, atol=1e-10,
        err_msg="Re-extraction of legacy flat with new code does not match legacy HDF5"
    )


@pytest.mark.preprocessed_data
def test_dragons_flat_self_consistency(
    preprocessed_files_path
):
    """
    Test 7: Box-extract the DRAGONS flat 2D image using new pipeline code
    and compare to what boxExtraction produces (BOX_REDUCED_FLAT_2).

    This validates that the pipeline's own extraction is self-consistent:
    manually calling _extract_single_stripe + _box_extract_single_stripe
    on the processed flat should match BOX_REDUCED_FLAT_2.
    """
    arm = 'BLUE'
    calib_root = preprocessed_files_path / 'calibrations'
    flat_path = str(calib_root / _CALIB_FILES[arm]['flat'])
    dark_path = str(calib_root / _CALIB_FILES[arm]['dark'])

    raw_file = preprocessed_files_path / (ETALON + '.fits')
    adinput = [astrodata.open(raw_file)]

    p = MaroonXSpectrum(adinput)
    p.prepare()
    p.checkArm()
    p.addDQ()
    p.overscanCorrect()
    p.correctImageOrientation()
    p.addVAR(read_noise=True, poisson_noise=True)

    p.extractStripes(flat=flat_path, dark=dark_path)

    # Grab the mask and flat stripe BEFORE boxExtraction (which may not
    # preserve STRIPES_MASKS on the output)
    ad_pre = p.streams['main'][0]
    fiber_key = f'fiber_{FIBER}'
    science_mask = ad_pre[0].STRIPES_MASKS[fiber_key][ORDER]
    flat_stripe_from_pipeline = ad_pre[0].F_STRIPES[fiber_key][ORDER]

    adout = p.boxExtraction()
    ad = adout[0]

    # Get pipeline's BOX_REDUCED_FLAT_2
    new_orders = ad[0].REDUCED_ORDERS_FIBER_2.astype(int)
    order_idx = list(new_orders).index(int(ORDER))
    pipeline_flat_box = ad[0].BOX_REDUCED_FLAT_2[order_idx]

    # Manually box-extract using the same stripe and mask the pipeline used
    manual_flat_box = _box_extract_single_stripe(
        flat_stripe_from_pipeline, science_mask
    )
    # Apply the same zero→NaN conversion that boxExtraction now does
    manual_flat_box[manual_flat_box == 0] = np.nan

    # Compare only valid (non-NaN) pixels
    both_valid = ~np.isnan(pipeline_flat_box) & ~np.isnan(manual_flat_box)
    diff = np.abs(pipeline_flat_box[both_valid] - manual_flat_box[both_valid])

    print(f"\n--- Test 7: DRAGONS flat self-consistency ---")
    print(f"  Pipeline BOX_REDUCED_FLAT_2 shape: {pipeline_flat_box.shape}")
    print(f"  Manual extraction shape: {manual_flat_box.shape}")
    print(f"  NaN pixels: pipeline={np.isnan(pipeline_flat_box).sum()}, "
          f"manual={np.isnan(manual_flat_box).sum()}")
    print(f"  Valid pixels compared: {both_valid.sum()} / {len(pipeline_flat_box)}")
    print(f"  Max abs diff: {diff.max():.6e}" if len(diff) > 0 else "  No valid pixels")
    print(f"  Pixels differing (>1e-10): {np.sum(diff > 1e-10)}")

    if np.sum(diff > 1e-10) == 0:
        print("\n  -> CONCLUSION: Pipeline extraction is self-consistent.")
        print("     Manual extraction matches BOX_REDUCED_FLAT_2 exactly.")

    np.testing.assert_allclose(
        pipeline_flat_box[both_valid], manual_flat_box[both_valid],
        rtol=0, atol=1e-10,
        err_msg="Pipeline BOX_REDUCED_FLAT_2 does not match manual extraction"
    )


@pytest.mark.preprocessed_data
def test_dtype_explains_box_extraction_diff(
    path_to_legacy_flats, preprocessed_files_path
):
    """
    Test 8: Test whether dtype difference explains the flat box extraction
    discrepancy.

    Strategy: load the DRAGONS flat, cast it to float32, extract a stripe,
    box-sum it, and compare to the legacy HDF5 box extraction.

    Finding: The f64 and f32 diffs are identical for all orders (~5-38 counts),
    proving that the discrepancy is NOT caused by dtype but by genuine
    differences in the 2D flat processing (overscan, straylight, combination).
    """
    from astropy.io import fits

    arm = 'BLUE'
    calib_root = preprocessed_files_path / 'calibrations'
    dragons_flat_path = str(calib_root / _CALIB_FILES[arm]['flat'])
    legacy_hdf_path = str(
        path_to_legacy_flats / LEGACY_FFFFF_HDF
    )

    # Load DRAGONS flat and cast to float32 (matching legacy precision)
    dragons_ad = astrodata.open(dragons_flat_path)
    dragons_data_f32 = dragons_ad[0].data.astype(np.float32)

    # Get polynomial coefficients
    p_id_table = dragons_ad[0].STRIPES_ID
    p_id = {}
    for i, f in enumerate(dragons_ad[0].STRIPES_FIBERS):
        p_id[f'fiber_{f}'] = dict(
            (colname, p_id_table[6 * i : 6 * (i + 1)][colname].data)
            for colname in p_id_table.colnames
        )

    # Load legacy config BPM
    config_path = (
        path_to_legacy_flats.parent.parent.parent
        / 'Maroonx_configfiles'
        / '202411xx'
        / 'config_b.hdf'
    )
    import h5py as h5
    with h5.File(str(config_path), 'r') as cf:
        legacy_bpm_raw = np.array(cf['bad_pixel_map'])
        valid = np.array(cf['valid'])

    legacy_bpm = legacy_bpm_raw[valid == 1].reshape(
        (np.sum(valid[:, 2000]), np.sum(valid[2000]))
    )
    legacy_bpm = np.flip(legacy_bpm, (0, 1))

    slit_height = 10
    fiber_key = f'fiber_{FIBER}'

    # Orders 113 and 116 have a known NaN discrepancy: the legacy flat has
    # 1 NaN pixel in these orders while the new pipeline has 0.0. This causes
    # an integer-level box-sum difference (~27, ~25 counts) that is unrelated
    # to the dtype issue being tested here. We track them separately.
    KNOWN_NAN_ORDERS = {'113', '116'}

    print(f"\n--- Test 8: dtype explains box extraction discrepancy ---")
    print(f"{'Order':>6} | {'f64 diff':>10} | {'f32 diff':>10} | {'f32 match':>10} | {'note':>8}")
    print("-" * 65)

    n_orders_f32_match = 0
    n_orders_total = 0
    nan_order_diffs = []

    for order_key, poly_coeffs in p_id[fiber_key].items():
        n_orders_total += 1

        # Extract at float64 (what DRAGONS does)
        stripe_f64 = _extract_single_stripe(
            dragons_ad[0].data, poly_coeffs, slit_height
        )
        bpm_stripe = _extract_single_stripe(
            legacy_bpm.astype(int), poly_coeffs, slit_height
        )
        box_f64 = _box_extract_single_stripe(stripe_f64, bpm_stripe)

        # Extract at float32 (matching legacy)
        stripe_f32 = _extract_single_stripe(
            dragons_data_f32, poly_coeffs, slit_height
        )
        box_f32 = _box_extract_single_stripe(stripe_f32, bpm_stripe)

        # Legacy box extraction from HDF5
        with tables.open_file(legacy_hdf_path, 'r') as h5flat:
            legacy_box = np.array(
                h5flat.get_node(f'/box_extraction/fiber_{FIBER}/{order_key}')
            )

        both_valid = ~np.isnan(box_f64) & ~np.isnan(legacy_box)

        diff_f64 = np.max(np.abs(
            box_f64[both_valid] - legacy_box[both_valid]
        )) if both_valid.sum() > 0 else 0.0

        diff_f32 = np.max(np.abs(
            box_f32[both_valid] - legacy_box[both_valid]
        )) if both_valid.sum() > 0 else 0.0

        is_nan_order = order_key in KNOWN_NAN_ORDERS
        f32_match = diff_f32 < 1e-10

        if f32_match:
            n_orders_f32_match += 1
        elif is_nan_order:
            nan_order_diffs.append((order_key, diff_f32))

        note = 'NaN' if is_nan_order else ''
        print(
            f"{order_key:>6} | {diff_f64:>10.2e} | {diff_f32:>10.2e} | "
            f"{'YES' if f32_match else 'NO':>10} | {note:>8}"
        )

    n_non_nan = n_orders_total - len(KNOWN_NAN_ORDERS)
    print(f"\n  Orders matching at float32: {n_orders_f32_match} / {n_orders_total}")
    print(f"  Excluding known NaN orders ({', '.join(sorted(KNOWN_NAN_ORDERS))}): "
          f"{n_orders_f32_match} / {n_non_nan}")

    if nan_order_diffs:
        print(f"\n  NaN-affected orders (separate issue — legacy has NaN, new has 0.0):")
        for o, d in nan_order_diffs:
            print(f"    Order {o}: diff = {d:.2e}")

    if n_orders_f32_match >= n_non_nan:
        print("\n  -> CONCLUSION: All non-NaN orders match when extracted at float32.")
        print("     The flat box extraction discrepancy is 100% explained by")
        print("     the DRAGONS flat being stored as float64 vs legacy float32.")
        print("     Orders 113/116 differ due to a separate NaN handling issue.")

    # All non-NaN orders must match at float32
    assert n_orders_f32_match >= n_non_nan, (
        f"Only {n_orders_f32_match}/{n_non_nan} non-NaN orders match at float32. "
        f"There may be processing differences beyond dtype."
    )
