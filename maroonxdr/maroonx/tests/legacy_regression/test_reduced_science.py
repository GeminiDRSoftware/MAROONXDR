from pathlib import Path

import astrodata
import numpy as np
import pandas as pd
from gempy.adlibrary import dataselect
from gempy.utils import logutils

import maroonx_instruments  # noqa : important to load adclass tags
from maroonxdr.maroonx.primitives_maroonx_spectrum import MaroonXSpectrum

#from maroonxdr.maroonx.primitives_maroonx_echelle import create_synthetic_dark

# Set logger
logutils.config(file_name="test_reduced_science.log", mode="debug", stomp=True)
log = logutils.get_logger("test_reduced_science.log")
log.setLevel("DEBUG")

# =========================================================
# TESTS
# =========================================================


def test_exportBundle(legacy_reduced_path, preprocessed_files_path):

    # Read legacy HDF5 file for comparison (300s SOOOE_x exposure)
    legacy_file = legacy_reduced_path / "20241124T062858Z_SOOOE_x_0300.hd5"

    # Read processed FITS files from DRAGONS reduction
    # These should be the output from a complete reduction workflow
    files = sorted([str(f) for f in preprocessed_files_path.glob('20241124T062858Z*.fits')])
    selected_spect = dataselect.select_data(files,
        tags=['PROCESSED', 'SCI', '300s'], xtags=['BUNDLE'])

    # Separate arms into streams and re-bundle
    adinput = [astrodata.open(f) for f in selected_spect]
    p = MaroonXSpectrum(adinput)
    p.separateArmStreams()
    bundled_list = p.bundleArmStreams()

    assert len(bundled_list) == 1, 'Should produce one bundle'
    bundle_ad = bundled_list[0]

    # Verify bundle structure matches what would be needed for HDF5 export
    assert len(bundle_ad) == 2, 'Bundle should have 2 extensions'
    assert 'BLUE' in bundle_ad[0].tags, 'First extension should be BLUE arm'
    assert 'RED' in bundle_ad[1].tags, 'Second extension should be RED arm'

    # Get both arms from bundle
    blue_ad, red_ad = bundle_ad[:]

    # Check that both arms have extracted spectral data for expected fibers
    expected_fibers = [2, 3, 4, 6]

    for fiber in expected_fibers:
        # Check BLUE arm attributes
        assert hasattr(blue_ad, f'OPTIMAL_REDUCED_FIBER_{fiber}'), \
            f'BLUE arm missing OPTIMAL_REDUCED_FIBER_{fiber}'

        # Check RED arm attributes
        assert hasattr(red_ad, f'OPTIMAL_REDUCED_FIBER_{fiber}'), \
            f'RED arm missing OPTIMAL_REDUCED_FIBER_{fiber}'

    # Compare with legacy HDF5 structure
    # Legacy files use pandas HDFStore with tables: 'spec_blue' and 'spec_red'
    # Each table has MultiIndex (fiber, order) and columns:
    # ['box_extraction', 'wavelengths', 'optimal_extraction', 'optimal_var']
    with pd.HDFStore(legacy_file, 'r') as legacy_store:

        legacy_blue = legacy_store['spec_blue']
        legacy_red = legacy_store['spec_red']

        # Compare BLUE arm data
        blue_failures = 0
        for fiber in expected_fibers:
            # Get DRAGONS data for this fiber
            dragons_wls = getattr(blue_ad, f'WLS_SIMULTANEOUS_FIBER_{fiber}')
            dragons_flux = getattr(blue_ad, f'OPTIMAL_REDUCED_FIBER_{fiber}')
            dragons_orders = blue_ad.REDUCED_ORDERS_FIBER_3.astype(int)

            # Compare each order
            for idx, order in enumerate(dragons_orders):
                if (fiber, order) not in legacy_blue.index:
                    log.warning(f"Fiber {fiber}, order {order} not in legacy BLUE data")
                    continue

                legacy_row = legacy_blue.loc[(fiber, order)]
                legacy_wls = legacy_row['wavelengths']
                legacy_flux = legacy_row['optimal_extraction']

                # Compare wavelengths
                try:
                    np.testing.assert_allclose(
                        dragons_wls[idx, :], legacy_wls,
                        rtol=0, atol=1e-4,
                        err_msg=f"BLUE fiber {fiber} order {order} wavelength mismatch"
                    )
                except AssertionError:
                    log.error(f"BLUE fiber {fiber} order {order} wavelength: FAIL")
                    blue_failures += 1

                # Compare flux
                try:
                    np.testing.assert_allclose(
                        dragons_flux[idx, :], legacy_flux,
                        rtol=1e-3, atol=1e-4,
                        err_msg=f"BLUE fiber {fiber} order {order} flux mismatch"
                    )
                except AssertionError:
                    log.error(f"BLUE fiber {fiber} order {order} flux: FAIL")
                    blue_failures += 1

        # Compare RED arm data
        red_failures = 0
        for fiber in expected_fibers:
            # Get DRAGONS data for this fiber
            dragons_wls = getattr(red_ad, f'WLS_SIMULTANEOUS_FIBER_{fiber}')
            dragons_flux = getattr(red_ad, f'OPTIMAL_REDUCED_FIBER_{fiber}')
            dragons_orders = red_ad.REDUCED_ORDERS_FIBER_3.astype(int)

            # Compare each order
            for idx, order in enumerate(dragons_orders):
                if (fiber, order) not in legacy_red.index:
                    log.warning(f"Fiber {fiber}, order {order} not in legacy RED data")
                    continue

                legacy_row = legacy_red.loc[(fiber, order)]
                legacy_wls = legacy_row['wavelengths']
                legacy_flux = legacy_row['optimal_extraction']

                # Compare wavelengths
                try:
                    np.testing.assert_allclose(
                        dragons_wls[idx, :], legacy_wls,
                        rtol=0, atol=1e-4,
                        err_msg=f"RED fiber {fiber} order {order} wavelength mismatch"
                    )
                except AssertionError:
                    log.error(f"RED fiber {fiber} order {order} wavelength: FAIL")
                    red_failures += 1

                # Compare flux
                try:
                    np.testing.assert_allclose(
                        dragons_flux[idx, :], legacy_flux,
                        rtol=1e-3, atol=1e-4,
                        err_msg=f"RED fiber {fiber} order {order} flux mismatch"
                    )
                except AssertionError:
                    log.error(f"RED fiber {fiber} order {order} flux: FAIL")
                    red_failures += 1

        # Report results
        total_failures = blue_failures + red_failures
        assert total_failures == 0, f"Found {total_failures} mismatches between DRAGONS and legacy data"

