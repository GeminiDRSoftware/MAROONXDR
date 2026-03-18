import os
from pathlib import Path

import astrodata
import numpy as np
import pytest
from gempy.utils import logutils

import maroonx_instruments  # noqa : important to load adclass tags
from maroonxdr.maroonx.tests.conftest import assert_allclose_with_max_fails

from . import legacy_adapter

# Set logger
logutils.config(file_name="test_reduced_wavecal.log", mode="debug", stomp=True)
log = logutils.get_logger("test_reduced_wavecal.log")
log.setLevel("DEBUG")

# =========================================================
# TESTS
# =========================================================


@pytest.mark.slow
@pytest.mark.preprocessed_data
@pytest.mark.parametrize(
    'matching_filenames',
    [
        # (new DRAGONS file, legacy HDF5 file) — etalon frames
        ('20241124T030227Z_DEEEE_b_0030_wavecal.fits', '20241124T030227Z_DEEEE_b_0030.hdf'),
        ('20241124T030227Z_DEEEE_r_0004_wavecal.fits', '20241124T030227Z_DEEEE_r_0004.hdf'),
    ],
)
def test_etalon_box_extraction(path_to_inputs, path_to_legacy_wavecal, matching_filenames):

    file_name, legacy_file_name = matching_filenames

    legacy_file = path_to_legacy_wavecal / legacy_file_name
    file = os.path.join(path_to_inputs, file_name)

    ad = astrodata.open(file)

    legacy_box = legacy_adapter.load_dict_from_hdf5(str(legacy_file), 'box_extraction/')

    for fiber in range(1, 6):
        new_box = getattr(ad[0], f'BOX_REDUCED_FIBER_{fiber}')
        if new_box.size == 1:
            continue

        orders = ad[0].REDUCED_ORDERS_FIBER_3.astype(int)

        for idx, order in enumerate(orders):
            legacy_order_box = legacy_box[f'fiber_{fiber}'][f'{order}']

            np.testing.assert_allclose(
                new_box[idx, :], legacy_order_box,
                rtol=1e-3, atol=1e-4,
                err_msg=f'fiber {fiber} order {order} box_extraction mismatch')




@pytest.mark.slow
@pytest.mark.preprocessed_data
@pytest.mark.parametrize(
    'matching_filenames',
    [
        # (new DRAGONS file, legacy HDF5 file) — etalon frames
        ('20241124T030227Z_DEEEE_b_0030_wavecal.fits', '20241124T030227Z_DEEEE_b_0030.hdf'),
        ('20241124T030227Z_DEEEE_r_0004_wavecal.fits', '20241124T030227Z_DEEEE_r_0004.hdf'),
    ],
)
def test_etalon_wavelengths(path_to_inputs, path_to_legacy_wavecal, matching_filenames):

    file_name, legacy_file_name = matching_filenames

    legacy_file = path_to_legacy_wavecal / legacy_file_name
    file = os.path.join(path_to_inputs, file_name)

    ad = astrodata.open(file)

    legacy_wls = legacy_adapter.load_dict_from_hdf5(str(legacy_file), 'wavelengths/')

    for fiber in range(1, 6):
        new_wls = getattr(ad[0], f'WLS_DYNAMIC_FIBER_{fiber}')
        if new_wls.size == 1:
            continue

        orders = ad[0].REDUCED_ORDERS_FIBER_3.astype(int)

        for idx, order in enumerate(orders):
            legacy_order_wls = legacy_wls[f'fiber_{fiber}'][f'{order}']

            np.testing.assert_allclose(
                new_wls[idx, :], legacy_order_wls,
                rtol=1e-3, atol=1e-4,
                err_msg=f'fiber {fiber} order {order} wavelength mismatch')


@pytest.mark.slow
@pytest.mark.preprocessed_data
@pytest.mark.parametrize(
    'matching_filenames',
    [
        # (new DRAGONS file, legacy HDF5 file) — LFC frames
        ('20241124T030436Z_DLLLL_b_0005_wavecal.fits', '20241124T030436Z_DLLLL_b_0005.hdf'),
        ('20241124T030436Z_DLLLL_r_0004_wavecal.fits', '20241124T030436Z_DLLLL_r_0004.hdf'),
    ],
)
def test_lfc_box_extraction(path_to_inputs, path_to_legacy_wavecal, matching_filenames):

    file_name, legacy_file_name = matching_filenames

    legacy_file = path_to_legacy_wavecal / legacy_file_name
    file = os.path.join(path_to_inputs, file_name)

    ad = astrodata.open(file)

    legacy_box = legacy_adapter.load_dict_from_hdf5(str(legacy_file), 'box_extraction/')

    for fiber in range(1, 6):
        new_box = getattr(ad[0], f'BOX_REDUCED_FIBER_{fiber}')
        if new_box.size == 1:
            continue

        orders = ad[0].REDUCED_ORDERS_FIBER_3.astype(int)

        for idx, order in enumerate(orders):
            legacy_order_box = legacy_box[f'fiber_{fiber}'][f'{order}']

            np.testing.assert_allclose(
                new_box[idx, :], legacy_order_box,
                rtol=1e-3, atol=1e-4,
                err_msg=f'fiber {fiber} order {order} box_extraction mismatch')




@pytest.mark.slow
@pytest.mark.preprocessed_data
@pytest.mark.parametrize(
    'matching_filenames',
    [
        # (new DRAGONS file, legacy HDF5 file) — LFC frames
        ('20241124T030436Z_DLLLL_b_0005_wavecal.fits', '20241124T030436Z_DLLLL_b_0005.hdf'),
        ('20241124T030436Z_DLLLL_r_0004_wavecal.fits', '20241124T030436Z_DLLLL_r_0004.hdf'),
    ],
)
def test_lfc_wavelengths(path_to_inputs, path_to_legacy_wavecal, matching_filenames):

    file_name, legacy_file_name = matching_filenames

    legacy_file = path_to_legacy_wavecal / legacy_file_name
    file = os.path.join(path_to_inputs, file_name)

    ad = astrodata.open(file)

    legacy_wls = legacy_adapter.load_dict_from_hdf5(str(legacy_file), 'wavelengths/')

    for fiber in range(1, 6):
        new_wls = getattr(ad[0], f'WLS_DYNAMIC_FIBER_{fiber}')
        if new_wls.size == 1:
            continue

        orders = ad[0].REDUCED_ORDERS_FIBER_3.astype(int)

        for idx, order in enumerate(orders):
            legacy_order_wls = legacy_wls[f'fiber_{fiber}'][f'{order}']

            np.testing.assert_allclose(
                new_wls[idx, :], legacy_order_wls,
                rtol=1e-3, atol=1e-4,
                err_msg=f'fiber {fiber} order {order} wavelength mismatch')
