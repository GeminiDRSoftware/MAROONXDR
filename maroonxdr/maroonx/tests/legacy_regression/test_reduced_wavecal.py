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

    fail_counter = 0
    for fiber in range(1, 6):
        new_box = getattr(ad[0], f'BOX_REDUCED_FIBER_{fiber}')
        if new_box.size == 1:
            continue

        orders = getattr(ad[0], f'REDUCED_ORDERS_FIBER_{fiber}').astype(int)

        for idx, order in enumerate(orders):
            legacy_order_box = legacy_box[f'fiber_{fiber}'][f'{order}']

            try:
                assert_allclose_with_max_fails(
                    new_box[idx, :], legacy_order_box,
                    rtol=0, atol=1e-8, max_fails=0)
                log.fullinfo(f'fiber/order : {fiber}/{order} [OK]')
            except (AssertionError, pytest.xfail.Exception):
                fail_counter += 1
                log.fullinfo(f'fiber/order : {fiber}/{order} [FAIL]')

    assert fail_counter == 0, f"{fail_counter} failing f/o combinations"


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

    fail_counter = 0
    for fiber in range(1, 6):
        new_wls = getattr(ad[0], f'WLS_DYNAMIC_FIBER_{fiber}')
        if new_wls.size == 1:
            continue

        orders = getattr(ad[0], f'REDUCED_ORDERS_FIBER_{fiber}').astype(int)

        for idx, order in enumerate(orders):
            legacy_order_wls = legacy_wls[f'fiber_{fiber}'][f'{order}']

            try:
                assert_allclose_with_max_fails(
                    new_wls[idx, :], legacy_order_wls,
                    rtol=0, atol=1e-8, max_fails=0)
                log.fullinfo(f'fiber/order : {fiber}/{order} [OK]')
            except (AssertionError, pytest.xfail.Exception):
                fail_counter += 1
                log.fullinfo(f'fiber/order : {fiber}/{order} [FAIL]')

    assert fail_counter == 0, f"{fail_counter} failing f/o combinations"


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

    if not os.path.exists(file):
        pytest.skip(f'LFC wavecal input not found: {file_name}')

    ad = astrodata.open(file)

    legacy_box = legacy_adapter.load_dict_from_hdf5(str(legacy_file), 'box_extraction/')

    fail_counter = 0
    for fiber in range(1, 6):
        new_box = getattr(ad[0], f'BOX_REDUCED_FIBER_{fiber}')
        if new_box.size == 1:
            continue

        orders = getattr(ad[0], f'REDUCED_ORDERS_FIBER_{fiber}').astype(int)

        for idx, order in enumerate(orders):
            legacy_order_box = legacy_box[f'fiber_{fiber}'][f'{order}']

            try:
                assert_allclose_with_max_fails(
                    new_box[idx, :], legacy_order_box,
                    rtol=0, atol=1e-8, max_fails=0)
                log.fullinfo(f'fiber/order : {fiber}/{order} [OK]')
            except (AssertionError, pytest.xfail.Exception):
                fail_counter += 1
                log.fullinfo(f'fiber/order : {fiber}/{order} [FAIL]')

    assert fail_counter == 0, f"{fail_counter} failing f/o combinations"


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

    if not os.path.exists(file):
        pytest.skip(f'LFC wavecal input not found: {file_name}')

    ad = astrodata.open(file)

    legacy_wls = legacy_adapter.load_dict_from_hdf5(str(legacy_file), 'wavelengths/')

    fail_counter = 0
    for fiber in range(1, 6):
        new_wls = getattr(ad[0], f'WLS_DYNAMIC_FIBER_{fiber}')
        if new_wls.size == 1:
            continue

        orders = getattr(ad[0], f'REDUCED_ORDERS_FIBER_{fiber}').astype(int)

        for idx, order in enumerate(orders):
            legacy_order_wls = legacy_wls[f'fiber_{fiber}'][f'{order}']

            try:
                assert_allclose_with_max_fails(
                    new_wls[idx, :], legacy_order_wls,
                    rtol=0, atol=1e-8, max_fails=0)
                log.fullinfo(f'fiber/order : {fiber}/{order} [OK]')
            except (AssertionError, pytest.xfail.Exception):
                fail_counter += 1
                log.fullinfo(f'fiber/order : {fiber}/{order} [FAIL]')

    assert fail_counter == 0, f"{fail_counter} failing f/o combinations"
