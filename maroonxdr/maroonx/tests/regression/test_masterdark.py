import re
import numpy as np
import pytest
from pathlib import Path

from astropy.io import fits

# from test_utils import compare_npz_files, compare_fits_files, compare_functions


# =========================================================
# FILES TO COMPARE
# =========================================================

# Paths to masterdarks coefficient files
# COEFF_OUTPUT_PATH = Path("/home/martin/Documentos/Projects/MaroonX/maroonx_base/data2/MaroonX_spectra_reduced/Maroonx_masterframes/202411xx/darks/")
# COEFF_REFERENCE_PATH = Path("/home/martin/Documentos/Projects/MaroonX/resources/Martín_setup/darks/")

# Paths to real masterdarks files
OLD_FILES_PATH = Path("/home/martin/Documentos/Projects/MaroonX/maroonx_base/data2/MaroonX_spectra_reduced/Maroonx_masterframes/202411xx/darks")
NEW_FILES_PATH = Path("/home/martin/Documentos/Projects/MAROONXDR/calibrations/processed_dark")

# Paths to synth masterdarks files
# SYNTH_OUTPUT_PATH = Path("/home/martin/Documentos/Projects/MaroonX/maroonx_base/data2/MaroonX_spectra_reduced/Maroonx_masterframes/202411xx/darks/")
# SYNTH_REFERENCE_PATH = Path("/home/martin/Documentos/Projects/MaroonX/resources/Martín_setup/darks/202411xx_synth_darks/darks/")

# Lists of all file names produced by the reduction process
# COEFF_FILES = [
#     "masterdarks_202411xx_blue.npz",
#     "masterdarks_202411xx_red.npz",
#     "masterdarks_coeffs_202411xx_blue.npz",
#     "masterdarks_coeffs_202411xx_red.npz"
# ]

# REAL_MASTERDARK_FILES = [
#     "20241115T19_masterdark_mean_DDDDE_b_0060.fits",
#     "20241115T19_masterdark_mean_DDDDE_b_0120.fits",
#     "20241115T19_masterdark_mean_DDDDE_b_0300.fits",
#     "20241115T19_masterdark_mean_DDDDE_r_0060.fits",
#     "20241115T19_masterdark_mean_DDDDE_r_0120.fits",
#     "20241115T19_masterdark_mean_DDDDE_r_0300.fits",
#     "20241115T20_masterdark_mean_DDDDE_b_0600.fits",
#     "20241115T20_masterdark_mean_DDDDE_r_0600.fits",
#     "20241115T21_masterdark_mean_DDDDE_b_0900.fits",
#     "20241115T21_masterdark_mean_DDDDE_r_0900.fits",
#     "20241115T22_masterdark_mean_DDDDE_b_1200.fits",
#     "20241115T22_masterdark_mean_DDDDE_r_1200.fits",
#     "20241116T00_masterdark_mean_DDDDE_b_1800.fits",
#     "20241116T00_masterdark_mean_DDDDE_r_1800.fits",
# ]

# SYNTH_MASTERDARK_FILES = [
#     "202411xx_masterdark_mean_DDDDE_b_0120.fits",
#     "202411xx_masterdark_mean_DDDDE_b_0300.fits",
#     "202411xx_masterdark_mean_DDDDE_b_0600.fits",
#     "202411xx_masterdark_mean_DDDDE_b_0900.fits",
#     "202411xx_masterdark_mean_DDDDE_b_1800.fits",
#     "202411xx_masterdark_mean_DDDDE_r_0120.fits",
#     "202411xx_masterdark_mean_DDDDE_r_0300.fits",
#     "202411xx_masterdark_mean_DDDDE_r_0600.fits",
#     "202411xx_masterdark_mean_DDDDE_r_0900.fits",
#     "202411xx_masterdark_mean_DDDDE_r_1800.fits"    
# ]
# =========================================================
# HELPER FUNCTIONS
# =========================================================


# =========================================================
# TESTS
# =========================================================


# The following tests should be sorted in order of pipeline
# output order to make it easier to traceback errors

@pytest.mark.parametrize("old_filename", ["20241115T19_masterdark_mean_DDDDE_b_0300.fits"])
@pytest.mark.parametrize("new_filename", ["20241115T193254Z_DDDDE_b_0300_regressionDark.fits"])
def test_real_comparison(old_filename, new_filename):

    is_blue = "_b_" in old_filename

    old_file = OLD_FILES_PATH / old_filename
    new_file = NEW_FILES_PATH / new_filename

    with fits.open(old_file) as old_hdu, fits.open(new_file) as new_hdu:
        old_data = old_hdu[0].data
        new_data = new_hdu[1].data

        assert old_data.shape == new_data.shape, f"Shape mismatch: {old_data.shape} != {new_data.shape}"

        # Compare the data
        np.testing.assert_allclose(old_data, new_data, rtol=1e-5, atol=1e-5)