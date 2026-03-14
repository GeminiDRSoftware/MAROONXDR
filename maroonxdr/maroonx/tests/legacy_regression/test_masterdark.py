import os

import astrodata
import numpy as np
import pytest
from astropy.io import fits
from gempy.adlibrary import dataselect
from gempy.utils import logutils

import maroonx_instruments  # noqa : important to load adclass tags
from maroonxdr.maroonx.primitives_maroonx_2D import MAROONX


# Set logger
logutils.config(file_name="test_masterdarks.log", mode="debug", stomp=True)
log = logutils.get_logger("test_masterdarks.log")
log.setLevel("DEBUG")


@pytest.mark.slow
@pytest.mark.preprocessed_data
@pytest.mark.parametrize(
    'matching_filenames',
    [
        # Add more tuples of (new_file, legacy_file) as needed
        ('20241115T190028Z_DDDDE_r_0120_dark.fits', '20241115T19_masterdark_mean_DDDDE_r_0120.fits'),
        ('20241115T191909Z_DDDDE_b_0060_dark.fits','20241115T19_masterdark_mean_DDDDE_b_0060.fits'),
        ('20241115T210524Z_DDDDE_b_0900_dark.fits','20241115T21_masterdark_mean_DDDDE_b_0900.fits'),
        ('20241115T210524Z_DDDDE_r_0900_dark.fits','20241115T21_masterdark_mean_DDDDE_r_0900.fits'),
        ('20241116T001751Z_DDDDE_b_1800_dark.fits','20241116T00_masterdark_mean_DDDDE_b_1800.fits'),
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

        # Compare the data
        np.testing.assert_allclose(old_data, new_data, rtol=0, atol=1e-4,
            err_msg='Data mismatch')




@pytest.mark.slow()
@pytest.mark.preprocessed_data
@pytest.mark.parametrize(
    'matching_filenames',
    [
        # Add more tuples of (new_file, legacy_file) as needed
        ('20241115T190028Z_DDDDE_b_0120_darkCoefficients.fits', 'masterdarks_coeffs_202411xx_blue.npz'),
        ('20241115T190028Z_DDDDE_r_0120_darkCoefficients.fits', 'masterdarks_coeffs_202411xx_red.npz'),
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

    # Compare the data
    np.testing.assert_allclose(z0, legacy_z0, rtol=0, atol=2e-4, err_msg='Data mismatch')
    np.testing.assert_allclose(z1, legacy_z1, rtol=0, atol=2e-4, err_msg='Data mismatch')
    np.testing.assert_allclose(logexptime['logexptime'].value, legacy_coeffs['logexptime'])



@pytest.mark.slow()
@pytest.mark.preprocessed_data
@pytest.mark.parametrize(
    'matching_filenames',
    [
        # Add more tuples of (new_file, legacy_file) as needed
        ('20241124T041907Z_SOOOE_b_0300_synth_dark.fits', '202411xx_masterdark_mean_DDDDE_b_0300.fits'),
        ('20241124T075055Z_SOOOE_r_0900_synth_dark.fits', '202411xx_masterdark_mean_DDDDE_r_0900.fits'),
    ],
)
def test_synthetic_masterdark(path_to_inputs, path_to_legacy_darks, matching_filenames):
    """
    Test the creation of a synthetic master dark.
    """

    legacy_file = path_to_legacy_darks / matching_filenames[1]
    file = os.path.join(path_to_inputs, matching_filenames[0])

    ad_coeff = astrodata.open(file)

    with fits.open(legacy_file) as hdul:
        legacy_dark = hdul[0].data
        synthetic_dark = ad_coeff[0].data

        # check that the synthetic dark frame is within tolerance of the legacy dark frame
        np.testing.assert_allclose(synthetic_dark, legacy_dark, rtol=0, atol=1e-4)