import os

import astrodata
import numpy as np
import pytest
from astropy.io import fits

import maroonx_instruments  # noqa : important to load adclass tags


@pytest.mark.preprocessed_data
@pytest.mark.parametrize(
    'matching_filenames',
    [
        ('20250721T162823Z_DDDDE_r_0120_dark.fits', '20250721T16_masterdark_mean_DDDDE_r_0120.fits'),
        ('20250721T164703Z_DDDDE_b_0060_dark.fits', '20250721T16_masterdark_mean_DDDDE_b_0060.fits'),
        ('20250721T183319Z_DDDDE_b_0900_dark.fits', '20250721T18_masterdark_mean_DDDDE_b_0900.fits'),
        ('20250721T183319Z_DDDDE_r_0900_dark.fits', '20250721T18_masterdark_mean_DDDDE_r_0900.fits'),
        ('20250721T214546Z_DDDDE_b_1800_dark.fits', '20250721T21_masterdark_mean_DDDDE_b_1800.fits'),
    ],
)
def test_masterdark(path_to_inputs, path_to_legacy_darks, matching_filenames):

    file_name, legacy_file_name = matching_filenames

    legacy_file = path_to_legacy_darks / legacy_file_name
    file = os.path.join(path_to_inputs, file_name)

    ad = astrodata.open(file)

    with fits.open(legacy_file) as old_hdu:
        old_data = old_hdu[0].data
        new_data = ad[0].data

        assert old_data.shape == new_data.shape, f"Shape mismatch: {old_data.shape} != {new_data.shape}"

        np.testing.assert_allclose(old_data, new_data, rtol=0, atol=1e-15,
            err_msg='Data mismatch')


@pytest.mark.preprocessed_data
@pytest.mark.parametrize(
    'matching_filenames',
    [
        ('20250721T162823Z_DDDDE_b_0120_darkCoefficients.fits', 'masterdarks_coeffs_202507xx_blue.npz'),
        ('20250721T162823Z_DDDDE_r_0120_darkCoefficients.fits', 'masterdarks_coeffs_202507xx_red.npz'),
    ],
)
def test_dark_coeff(path_to_inputs, path_to_legacy_darks, matching_filenames):

    legacy_file = path_to_legacy_darks / matching_filenames[1]
    legacy_coeffs = np.load(legacy_file)

    ad = astrodata.open(os.path.join(path_to_inputs, matching_filenames[0]))

    z0 = ad[0].COEFF_Z0
    z1 = ad[0].COEFF_Z1
    logexptime = ad[0].LOGEXPTIME

    legacy_z0 = legacy_coeffs['z0']
    legacy_z1 = legacy_coeffs['z1']

    np.testing.assert_allclose(z0, legacy_z0, rtol=1e-5, atol=1e-15, err_msg='Data mismatch')
    np.testing.assert_allclose(z1, legacy_z1, rtol=1e-5, atol=1e-15, err_msg='Data mismatch')
    np.testing.assert_allclose(logexptime['logexptime'].value, legacy_coeffs['logexptime'])


@pytest.mark.preprocessed_data
@pytest.mark.parametrize(
    'matching_filenames',
    [
        ('20250717T144308Z_SOOOE_b_0300_synth_dark.fits', '202507xx_masterdark_mean_DDDDE_b_0300.fits'),
        ('20250717T144308Z_SOOOE_r_0300_synth_dark.fits', '202507xx_masterdark_mean_DDDDE_r_0300.fits'),
    ],
)
def test_synthetic_masterdark(path_to_inputs, path_to_legacy_darks, matching_filenames):

    legacy_file = path_to_legacy_darks / matching_filenames[1]
    file = os.path.join(path_to_inputs, matching_filenames[0])

    ad_coeff = astrodata.open(file)

    with fits.open(legacy_file) as hdul:
        legacy_dark = hdul[0].data
        synthetic_dark = ad_coeff[0].data

        np.testing.assert_allclose(synthetic_dark, legacy_dark, rtol=1e-5, atol=1e-15,
            err_msg='Data mismatch')
