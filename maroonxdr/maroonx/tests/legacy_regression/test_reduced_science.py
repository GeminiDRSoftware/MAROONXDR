import os

import astrodata
import numpy as np
import pytest
from gempy.utils import logutils

import maroonx_instruments  # noqa : important to load adclass tags
from maroonxdr.maroonx.tests.conftest import assert_allclose_with_max_fails

from . import legacy_adapter

logutils.config(file_name="test_reduced_science_v2.log", mode="debug", stomp=True)
log = logutils.get_logger("test_reduced_science_v2.log")
log.setLevel("DEBUG")


@pytest.mark.preprocessed_data
@pytest.mark.parametrize(
    'matching_filenames',
    [
        ('20250717T144308Z_SOOOE_b_0300_reduced.fits', '20250717T144308Z_SOOOE_b_0300.hdf'),
        ('20250717T144308Z_SOOOE_r_0300_reduced.fits', '20250717T144308Z_SOOOE_r_0300.hdf'),
    ],
)
def test_optimal_extraction(path_to_inputs, path_to_legacy_science, matching_filenames):

    file_name, legacy_file_name = matching_filenames

    legacy_file = path_to_legacy_science / legacy_file_name
    file = os.path.join(path_to_inputs, file_name)

    ad = astrodata.open(file)

    legacy_opt = legacy_adapter.load_dict_from_hdf5(str(legacy_file), 'optimal_extraction/')

    fail_counter = 0
    for fiber in [1, 2, 3, 4, 5]:
        new_opt = getattr(ad[0], f'OPTIMAL_REDUCED_FIBER_{fiber}', None)
        if new_opt is None or new_opt.size == 1:
            continue

        orders = getattr(ad[0], f'REDUCED_ORDERS_FIBER_{fiber}').astype(int)

        for idx, order in enumerate(orders):
            legacy_fiber = legacy_opt.get(f'fiber_{fiber}')
            if legacy_fiber is None:
                continue
            legacy_order_opt = legacy_fiber.get(f'{order}')
            if legacy_order_opt is None:
                continue

            try:
                assert_allclose_with_max_fails(
                    new_opt[idx, :], legacy_order_opt,
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
        ('20250717T144308Z_SOOOE_b_0300_reduced.fits', '20250717T144308Z_SOOOE_b_0300.hdf'),
        ('20250717T144308Z_SOOOE_r_0300_reduced.fits', '20250717T144308Z_SOOOE_r_0300.hdf'),
    ],
)
def test_optimal_var(path_to_inputs, path_to_legacy_science, matching_filenames):

    file_name, legacy_file_name = matching_filenames

    legacy_file = path_to_legacy_science / legacy_file_name
    file = os.path.join(path_to_inputs, file_name)

    ad = astrodata.open(file)

    legacy_var = legacy_adapter.load_dict_from_hdf5(str(legacy_file), 'optimal_var/')

    fail_counter = 0
    for fiber in [1, 2, 3, 4, 5]:
        new_var = getattr(ad[0], f'OPTIMAL_REDUCED_VAR_{fiber}', None)
        if new_var is None or new_var.size == 1:
            continue

        orders = getattr(ad[0], f'REDUCED_ORDERS_FIBER_{fiber}').astype(int)

        for idx, order in enumerate(orders):
            legacy_fiber = legacy_var.get(f'fiber_{fiber}')
            if legacy_fiber is None:
                continue
            legacy_order_var = legacy_fiber.get(f'{order}')
            if legacy_order_var is None:
                continue

            try:
                assert_allclose_with_max_fails(
                    new_var[idx, :], legacy_order_var,
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
        ('20250717T144308Z_SOOOE_b_0300_reduced.fits', '20250717T144308Z_SOOOE_b_0300.hdf'),
        ('20250717T144308Z_SOOOE_r_0300_reduced.fits', '20250717T144308Z_SOOOE_r_0300.hdf'),
    ],
)
def test_simultaneous_wavelengths(path_to_inputs, path_to_legacy_science, matching_filenames):

    file_name, legacy_file_name = matching_filenames

    legacy_file = path_to_legacy_science / legacy_file_name
    file = os.path.join(path_to_inputs, file_name)

    ad = astrodata.open(file)

    legacy_wls = legacy_adapter.load_dict_from_hdf5(str(legacy_file), 'wavelengths/')

    fail_counter = 0
    for fiber in [1, 2, 3, 4, 5]:
        new_wls = getattr(ad[0], f'WLS_SIMULTANEOUS_FIBER_{fiber}', None)
        if new_wls is None or new_wls.size == 1:
            continue

        orders = getattr(ad[0], f'REDUCED_ORDERS_FIBER_{fiber}').astype(int)

        for idx, order in enumerate(orders):
            legacy_fiber = legacy_wls.get(f'fiber_{fiber}')
            if legacy_fiber is None:
                continue
            legacy_order_wls = legacy_fiber.get(f'{order}')
            if legacy_order_wls is None:
                continue

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
        ('20250717T144308Z_SOOOE_b_0300_reduced.fits', '20250717T144308Z_SOOOE_b_0300.hdf'),
        ('20250717T144308Z_SOOOE_r_0300_reduced.fits', '20250717T144308Z_SOOOE_r_0300.hdf'),
    ],
)
def test_fiber_6(path_to_inputs, path_to_legacy_science, matching_filenames):

    file_name, legacy_file_name = matching_filenames

    legacy_file = path_to_legacy_science / legacy_file_name
    file = os.path.join(path_to_inputs, file_name)

    ad = astrodata.open(file)

    legacy_opt = legacy_adapter.load_dict_from_hdf5(str(legacy_file), 'optimal_extraction/')

    fail_counter = 0
    for fiber in [6]:
        new_opt = getattr(ad[0], f'OPTIMAL_REDUCED_FIBER_{fiber}', None)
        if new_opt is None or new_opt.size == 1:
            continue

        orders = getattr(ad[0], f'REDUCED_ORDERS_FIBER_{fiber}').astype(int)

        for idx, order in enumerate(orders):
            legacy_fiber = legacy_opt.get(f'fiber_{fiber}')
            if legacy_fiber is None:
                continue
            legacy_order_opt = legacy_fiber.get(f'{order}')
            if legacy_order_opt is None:
                continue

            try:
                # the tolerance 
                assert_allclose_with_max_fails(
                    new_opt[idx, :], legacy_order_opt,
                    rtol=1e-7, atol=0, max_fails=0)
                log.fullinfo(f'fiber/order : {fiber}/{order} [OK]')
            except (AssertionError, pytest.xfail.Exception):
                fail_counter += 1
                log.fullinfo(f'fiber/order : {fiber}/{order} [FAIL]')

    assert fail_counter == 0, f"{fail_counter} failing f/o combinations"