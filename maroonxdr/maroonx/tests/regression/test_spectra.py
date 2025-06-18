

from copy import deepcopy
import numpy as np
import pandas as pd
import pytest
from pathlib import Path

import astrodata
from astropy.io import fits
import h5py

from gempy.adlibrary import dataselect

import maroonx_instruments  # noqa : important to load adclass tags
from maroonxdr.maroonx.primitives_maroonx_spectrum import MaroonXSpectrum
from maroonxdr.maroonx.primitives_maroonx_2D import MAROONX


# =========================================================
# FILES TO COMPARE
# =========================================================

OLD_FILES_PATH = Path("/home/martin/Projects/MaroonX/legacy/maroonx_base/data2/MaroonX_spectra_reduced/20241124")
NEW_FILES_PATH = Path("/home/martin/Projects/MaroonX/MAROONXDR/science_dir")

SCIENCE_DIR = Path('/home/martin/Projects/MaroonX/MAROONXDR/science_dir')

# =========================================================
# TESTS
# =========================================================

# @pytest.mark.parametrize("old_filename", ["20241115T19_masterdark_mean_DDDDE_b_0300.fits"])
# @pytest.mark.parametrize("new_filename", ["20241115T193254Z_DDDDE_b_0300_regressionDark.fits"])
# def test_boxExtraction():
#     pass


@pytest.mark.parametrize("old_filename", ["20241115T19_masterdark_mean_DDDDE_b_0300.fits"])
@pytest.mark.parametrize("new_filename", ["20241115T193254Z_DDDDE_b_0300_regressionDark.fits"])
def test_getPeaks(old_filename, new_filename):

    old_file = OLD_FILES_PATH / old_filename
    new_file = NEW_FILES_PATH / new_filename

    # Load old peak data. columns are lowercase
    old_peak_data = pd.read_hdf(old_file, '/etalon_peak_parameters/peaks')

    # Load new peak data. columns are uppercase
    new_peak_data = astrodata.open(new_file)[0].PEAKS.to_pandas()

    # Test shapes are equal
    assert old_peak_data.shape == new_peak_data.shape

    # Test the the column FIBER has the same value counts
    assert new_peak_data["FIBER"].value_counts().equals(old_peak_data["fiber"].value_counts())