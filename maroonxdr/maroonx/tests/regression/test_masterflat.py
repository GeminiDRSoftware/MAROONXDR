from copy import deepcopy
import numpy as np
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

OLD_FILES_PATH = Path("/home/martin/Documentos/Projects/MaroonX/maroonx_base/data2/MaroonX_spectra_reduced/Maroonx_masterframes/202411xx/flats")
NEW_FILES_PATH = Path("/home/martin/Documentos/Projects/MAROONXDR/calibrations/processed_flat")

SCIENCE_DIR = Path('/home/martin/Documentos/Projects/MAROONXDR/science_dir')

# =========================================================
# TESTS
# =========================================================
# The following tests should be sorted in order of pipeline
# output order to make it easier to traceback errors

@pytest.mark.parametrize("old_filename,new_filename", [
    ("20241114T19_masterflat_backgroundsubtracted_FFFFF_b_0007.fits", "20241114T190714Z_DDDDF_b_0007_DFFFF_flat.fits"),
    ("20241114T19_masterflat_backgroundsubtracted_FFFFF_r_0002.fits", "20241114T190714Z_DDDDF_r_0002_DFFFF_flat.fits"),
])
def test_real_comparison(old_filename, new_filename):

    old_file = OLD_FILES_PATH / old_filename
    new_file = NEW_FILES_PATH / new_filename

    with fits.open(old_file) as old_hdu, fits.open(new_file) as new_hdu:
        old_data = old_hdu[0].data
        new_data = new_hdu[1].data

        assert old_data.shape == new_data.shape, f"Shape mismatch: {old_data.shape} != {new_data.shape}"


        # Compare the data
        np.testing.assert_allclose(old_data, new_data, rtol=1e-5, atol=1e-5)


@pytest.mark.parametrize("arm", ["BLUE", "RED"])
def test_stackFlats(arm):

    if arm == "BLUE":
        old_DFFFD = OLD_FILES_PATH / "20241114T18_masterflat_DFFFD_b_0008.fits"
        old_DDDDF = OLD_FILES_PATH / "20241114T19_masterflat_DDDDF_b_0007.fits"
    else:
        old_DFFFD = OLD_FILES_PATH / "20241114T18_masterflat_DFFFD_r_0002.fits"
        old_DDDDF = OLD_FILES_PATH / "20241114T19_masterflat_DDDDF_r_0002.fits"

    # get all flat files
    raw_files = sorted([str(f) for f in SCIENCE_DIR.glob('*.fits')])
    selected_spect = dataselect.select_data(raw_files, tags=['RAW', 'FLAT', arm])

    # read files and instantiate the primitive class
    adinput = [astrodata.open(f) for f in selected_spect]
    p = MAROONX(adinput)

    p.subtractOverscan()
    # p.trimOverscan()  # old files are not trimmed at this point
    # p.correctImageOrientation()  # blue files are not flipped at this point
    p.separateFlatStreams()  # creates 'DFFFD_flats' stream and leaves FDDDF flats in main stream

    # stack each stream separately
    p.stackFlats(stream='main', scale_mode='mean_frame', suffix='_DDDDF_flat')
    p.stackFlats(stream='DFFFD_flats', scale_mode='mean_frame', suffix='_DFFFD_flat')

    ad_DDDDF = p.streams['main'][0]
    ad_DFFFD = p.streams['DFFFD_flats'][0]

    # Compare with old masterflats
    with fits.open(old_DDDDF) as hdu_DDDDF, fits.open(old_DFFFD) as hdu_DFFFD:

        np.testing.assert_allclose(ad_DDDDF[0].data, hdu_DDDDF[0].data,
            rtol=1e-5, atol=1e-5)
        np.testing.assert_allclose(ad_DFFFD[0].data, hdu_DFFFD[0].data,
            rtol=1e-5, atol=1e-5)


@pytest.mark.parametrize("arm", ["BLUE", "RED"])
def test_identifyStripes(arm):

    if arm == "BLUE":
        old_DFFFD = OLD_FILES_PATH / "20241114T18_masterflat_DFFFD_b_0008.hdf"
        old_DDDDF = OLD_FILES_PATH / "20241114T19_masterflat_DDDDF_b_0007.hdf"
    else:
        old_DFFFD = OLD_FILES_PATH / "20241114T18_masterflat_DFFFD_r_0002.hdf"
        old_DDDDF = OLD_FILES_PATH / "20241114T19_masterflat_DDDDF_r_0002.hdf"

    # get all flat files
    raw_files = sorted([str(f) for f in SCIENCE_DIR.glob('*.fits')])
    selected_spect = dataselect.select_data(raw_files, tags=['RAW', 'FLAT', arm])

    # read files and instantiate the primitive class
    adinput = [astrodata.open(f) for f in selected_spect]
    p = MAROONX(adinput)

    p.subtractOverscan()
    p.trimOverscan()  # old files are not trimmed at this point
    p.correctImageOrientation()  # blue files are not flipped at this point
    p.separateFlatStreams()  # creates 'DFFFD_flats' stream and leaves FDDDF flats in main stream
    
    p.stackFlats(stream='main', scale_mode='mean_frame', suffix='_DDDDF_flat')
    p.stackFlats(stream='DFFFD_flats', scale_mode='mean_frame', suffix='_DFFFD_flat')

    p.findStripes()  # define stripe info to ultimately remove stray light in each stream
    p.findStripes(stream='DFFFD_flats')

    p.identifyStripes(selected_fibers='0,0,0,0,5') # identify stripes based on MX architecture files
    p.identifyStripes(stream='DFFFD_flats',selected_fibers='0,2,3,4,0')

    p_id_DDDDF = p.streams['main'][0][0].STRIPES_ID
    p_id_DFFFD = p.streams['DFFFD_flats'][0][0].STRIPES_ID

    old_p_id_DDDDF = load_dict_from_hdf5(str(old_DDDDF), 'extraction_parameters/')
    old_p_id_DFFFD = load_dict_from_hdf5(str(old_DFFFD), 'extraction_parameters/')

    for fiber in old_p_id_DDDDF.keys():
        for order in old_p_id_DDDDF[fiber].keys():
            np.testing.assert_allclose(p_id_DDDDF[fiber][order], old_p_id_DDDDF[fiber][order],
                rtol=1e-5, atol=1e-5)
            
    for fiber in old_p_id_DFFFD.keys():
        for order in old_p_id_DFFFD[fiber].keys():
            np.testing.assert_allclose(p_id_DFFFD[fiber][order], old_p_id_DFFFD[fiber][order],
                rtol=1e-5, atol=1e-5)



# =====================================================
# HDF5 helper functions
# =====================================================

def load_dict_from_hdf5(filename, path):
    if not path.endswith("/"):
        path += "/"
    if isinstance(filename, str):
        with h5py.File(filename, 'r', libver='latest') as h5file:
            return recursively_load_dict_contents_from_group(h5file, path)
    else:
        return recursively_load_dict_contents_from_group(filename, path)

def recursively_load_dict_contents_from_group(h5file, path):
    ans = {}
    for key, item in h5file[path].items():
        if isinstance(item, h5py._hl.dataset.Dataset):
            ans[key] = item[()]
        elif isinstance(item, h5py._hl.group.Group):
            ans[key] = recursively_load_dict_contents_from_group(h5file, path + key + '/')
    return ans