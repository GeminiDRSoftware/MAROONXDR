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
# TESTS
# =========================================================


def test_stackFlats(arm, legacy_flats_path):
    """
    Needs legacy output of step 1, create master flats: 
    reduce/recipes/make_master_flats.py
    """

    if arm == "BLUE":
        old_DFFFD = legacy_flats_path / "20241114T18_masterflat_DFFFD_b_0008.fits"
        old_DDDDF = legacy_flats_path / "20241114T19_masterflat_DDDDF_b_0007.fits"
    else:
        old_DFFFD = legacy_flats_path / "20241114T18_masterflat_DFFFD_r_0002.fits"
        old_DDDDF = legacy_flats_path / "20241114T19_masterflat_DDDDF_r_0002.fits"

    # get all flat files
    raw_files = sorted([str(f) for f in Path().glob('*.fits')])
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

def test_identifyStripes(arm, legacy_flats_path):
    """
    Needs legacy output of step 1, create master flats: 
    reduce/extraction.py
    """
    if arm == "BLUE":
        old_DFFFD = legacy_flats_path / "20241114T18_masterflat_DFFFD_b_0008.hdf"
        old_DDDDF = legacy_flats_path / "20241114T19_masterflat_DDDDF_b_0007.hdf"
    else:
        old_DFFFD = legacy_flats_path / "20241114T18_masterflat_DFFFD_r_0002.hdf"
        old_DDDDF = legacy_flats_path / "20241114T19_masterflat_DDDDF_r_0002.hdf"

    # get all flat files
    raw_files = sorted([str(f) for f in Path().glob('*.fits')])
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

def test_identifyStripes_legacyOrder(arm, legacy_flats_path):
    """
    Needs legacy output of step 1, create master flats: 
    reduce/extraction.py
    """
    if arm == "BLUE":
        old_DFFFD = legacy_flats_path / "20241114T18_masterflat_DFFFD_b_0008.hdf"
        old_DDDDF = legacy_flats_path / "20241114T19_masterflat_DDDDF_b_0007.hdf"
    else:
        old_DFFFD = legacy_flats_path / "20241114T18_masterflat_DFFFD_r_0002.hdf"
        old_DDDDF = legacy_flats_path / "20241114T19_masterflat_DDDDF_r_0002.hdf"

    # get all flat files
    raw_files = sorted([str(f) for f in Path().glob('*.fits')])
    selected_spect = dataselect.select_data(raw_files, tags=['RAW', 'FLAT', arm])

    # read files and instantiate the primitive class
    adinput = [astrodata.open(f) for f in selected_spect]
    p = MAROONX(adinput)

    p.prepare()
    p.checkArm()
    p.checkND()
    p.addDQ()
    p.addVAR(read_noise=True, poisson_noise=True)

    p.subtractOverscan()

    p.separateFlatStreams()  # creates 'DFFFD_flats' stream and leaves FDDDF flats in main stream

    p.stackFlats(stream='main', scale_mode='mean_frame', suffix='_DDDDF_flat')
    p.stackFlats(stream='DFFFD_flats', scale_mode='mean_frame', suffix='_DFFFD_flat')
    
    # Subtract overscan is run again in backgroundfit.py on the stacked flats
    p.subtractOverscan(stream='main')
    p.subtractOverscan(stream='DFFFD_flats')

    p.trimOverscan(stream='main')
    p.trimOverscan(stream='DFFFD_flats')

    p.correctImageOrientation(stream='main')
    p.correctImageOrientation(stream='DFFFD_flats')
    # ================================================
    
    p.findStripes(stream='main')  # define stripe info to ultimately remove stray light in each stream
    p.findStripes(stream='DFFFD_flats')

    p.identifyStripes(stream='main', selected_fibers='0,0,0,0,5') # identify stripes based on MX architecture files
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

#@pytest.mark.xfail(reason="Photutils Background2D changed from 1.02 to 2.02")
def test_removeStraylight(arm, legacy_flats_path):
    # This test does not subtract overscan twice as test_removeStraylight_legacyOrder does.

    if arm == "BLUE":
        old_DFFFD = legacy_flats_path / "20241114T18_masterflat_backgroundsubtracted_DFFFD_b_0008.fits"
        old_DDDDF = legacy_flats_path / "20241114T19_masterflat_backgroundsubtracted_DDDDF_b_0007.fits"
    else:
        old_DFFFD = legacy_flats_path / "20241114T18_masterflat_backgroundsubtracted_DFFFD_r_0002.fits"
        old_DDDDF = legacy_flats_path / "20241114T19_masterflat_backgroundsubtracted_DDDDF_r_0002.fits"

    # get all flat files
    raw_files = sorted([str(f) for f in Path().glob('*.fits')])
    selected_spect = dataselect.select_data(raw_files, tags=['RAW', 'FLAT', arm])

    # read files and instantiate the primitive class
    adinput = [astrodata.open(f) for f in selected_spect]
    p = MAROONX(adinput)

    p.prepare()
    p.checkArm()
    p.checkND()
    p.addDQ()
    p.addVAR(read_noise=True, poisson_noise=True)
    p.subtractOverscan()
    p.trimOverscan()  # old files are not trimmed at this point
    p.correctImageOrientation()  # blue files are not flipped at this point
    p.separateFlatStreams()  # creates 'DFFFD_flats' stream and leaves FDDDF flats in main stream
    
    p.stackFlats(stream='main', scale_mode='mean_frame', suffix='_flat')
    p.stackFlats(stream='DFFFD_flats', scale_mode='mean_frame', suffix='_flat')

    p.findStripes()  # define stripe info to ultimately remove stray light in each stream
    p.findStripes(stream='DFFFD_flats')

    p.identifyStripes(selected_fibers='0,0,0,0,5') # identify stripes based on MX architecture files
    p.identifyStripes(stream='DFFFD_flats',selected_fibers='0,2,3,4,0')

    p.defineFlatStripes()  # defines pixel inclusion for each flat region based on stripe ids
    p.defineFlatStripes(stream='DFFFD_flats')

    p.removeStrayLight(filter_size=19, box_size=20, snapshot=False)  # remove straylight from frame (this is why 2 partial illumination flat sets are necessary)
    p.removeStrayLight(stream='DFFFD_flats', filter_size=19, box_size=20, snapshot=False)

    ad_DDDDF = p.streams['main'][0]
    ad_DFFFD = p.streams['DFFFD_flats'][0]

    # Compare with old masterflats
    with fits.open(old_DDDDF) as hdu_DDDDF, fits.open(old_DFFFD) as hdu_DFFFD:

        np.testing.assert_allclose(ad_DDDDF[0].data, hdu_DDDDF[0].data,
            rtol=1e-5, atol=1e-5)
        np.testing.assert_allclose(ad_DFFFD[0].data, hdu_DFFFD[0].data,
            rtol=1e-5, atol=1e-5)

#@pytest.mark.xfail(reason="Photutils Background2D changed from 1.02 to 2.02")
def test_removeStraylight_legacyOrder(arm, legacy_flats_path):

    if arm == "BLUE":
        old_DFFFD = legacy_flats_path / "20241114T18_masterflat_backgroundsubtracted_DFFFD_b_0008.fits"
        old_DDDDF = legacy_flats_path / "20241114T19_masterflat_backgroundsubtracted_DDDDF_b_0007.fits"
    else:
        old_DFFFD = legacy_flats_path / "20241114T18_masterflat_backgroundsubtracted_DFFFD_r_0002.fits"
        old_DDDDF = legacy_flats_path / "20241114T19_masterflat_backgroundsubtracted_DDDDF_r_0002.fits"

    # get all flat files
    raw_files = sorted([str(f) for f in Path().glob('*.fits')])
    selected_spect = dataselect.select_data(raw_files, tags=['RAW', 'FLAT', arm])

    # read files and instantiate the primitive class
    adinput = [astrodata.open(f) for f in selected_spect]
    p = MAROONX(adinput)

    p.prepare()
    p.checkArm()
    p.checkND()
    p.addDQ()
    p.addVAR(read_noise=True, poisson_noise=True)
    p.subtractOverscan()

    p.separateFlatStreams()  # creates 'DFFFD_flats' stream and leaves FDDDF flats in main stream

    p.stackFlats(stream='main', scale_mode='mean_frame', suffix='_flat')
    p.stackFlats(stream='DFFFD_flats', scale_mode='mean_frame', suffix='_flat')
    
    # Subtract overscan is run again in backgroundfit.py on the stacked flats
    p.subtractOverscan(stream='main')
    p.subtractOverscan(stream='DFFFD_flats')

    p.trimOverscan(stream='main')
    p.trimOverscan(stream='DFFFD_flats')

    p.correctImageOrientation(stream='main')
    p.correctImageOrientation(stream='DFFFD_flats')
    # ================================================

    p.findStripes()  # define stripe info to ultimately remove stray light in each stream
    p.findStripes(stream='DFFFD_flats')

    p.identifyStripes(selected_fibers='0,0,0,0,5') # identify stripes based on MX architecture files
    p.identifyStripes(stream='DFFFD_flats',selected_fibers='0,2,3,4,0')

    p.defineFlatStripes()  # defines pixel inclusion for each flat region based on stripe ids
    p.defineFlatStripes(stream='DFFFD_flats')

    p.removeStrayLight(filter_size=19, box_size=20, snapshot=False)  # remove straylight from frame (this is why 2 partial illumination flat sets are necessary)
    p.removeStrayLight(stream='DFFFD_flats', filter_size=19, box_size=20, snapshot=False)

    ad_DDDDF = p.streams['main'][0]
    ad_DFFFD = p.streams['DFFFD_flats'][0]

    # Compare with old masterflats
    with fits.open(old_DDDDF) as hdu_DDDDF, fits.open(old_DFFFD) as hdu_DFFFD:

        np.testing.assert_allclose(ad_DDDDF[0].data, hdu_DDDDF[0].data,
            rtol=1e-5, atol=1e-5)
        np.testing.assert_allclose(ad_DFFFD[0].data, hdu_DFFFD[0].data,
            rtol=1e-5, atol=1e-5)

def test_combinedFlat(arm, legacy_flats_path):

    if arm == "BLUE":
        old_DFFFF = legacy_flats_path / "20241114T19_masterflat_backgroundsubtracted_FFFFF_b_0007.hdf"
    else:
        old_DFFFF = legacy_flats_path / "20241114T19_masterflat_backgroundsubtracted_FFFFF_r_0002.hdf"

    # get all flat files
    raw_files = sorted([str(f) for f in Path().glob('*.fits')])
    selected_spect = dataselect.select_data(raw_files, tags=['RAW', 'FLAT', arm])

    # read files and instantiate the primitive class
    adinput = [astrodata.open(f) for f in selected_spect]
    p = MAROONX(adinput)

    p.prepare()
    p.checkArm()
    p.checkND()
    p.addDQ()
    p.addVAR(read_noise=True, poisson_noise=True)
    p.subtractOverscan()
    p.separateFlatStreams()  # creates 'DFFFD_flats' stream and leaves FDDDF flats in main stream

    p.stackFlats(stream='main', scale_mode='mean_frame', suffix='_flat')
    p.stackFlats(stream='DFFFD_flats', scale_mode='mean_frame', suffix='_flat')
    
    # Subtract overscan is run again in backgroundfit.py on the stacked flats
    p.subtractOverscan(stream='main')
    p.subtractOverscan(stream='DFFFD_flats')
    
    p.trimOverscan(stream='main')
    p.trimOverscan(stream='DFFFD_flats')
    
    p.correctImageOrientation(stream='main')
    p.correctImageOrientation(stream='DFFFD_flats')
    
    p.findStripes()  # define stripe info to ultimately remove stray light in each stream
    p.findStripes(stream='DFFFD_flats')
    
    p.identifyStripes(selected_fibers='0,0,0,0,5') # identify stripes based on MX architecture files
    p.identifyStripes(stream='DFFFD_flats',selected_fibers='0,2,3,4,0')
    
    p.defineFlatStripes()  # defines pixel inclusion for each flat region based on stripe ids
    p.defineFlatStripes(stream='DFFFD_flats')
    
    p.removeStrayLight(filter_size=19, box_size=20, snapshot=False)  # remove straylight from frame (this is why 2 partial illumination flat sets are necessary)
    p.removeStrayLight(stream='DFFFD_flats', filter_size=19, box_size=20, snapshot=False)
    # =========================================================

    p.combineFlatStreams(stream='main', stream_2='DFFFD_flats')  # combine straylight-removed images
    p.clearStream(stream='DFFFD_flats') # remove second stream
    
    p.findStripes()  # re-run find/identify/define routine on combined frame
    p.identifyStripes(selected_fibers='0,2,3,4,5')
    # p.defineFlatStripes(extract=True)

    p_id_DFFFF = p.streams['main'][0][0].STRIPES_ID

    old_p_id_DFFFF = load_dict_from_hdf5(str(old_DFFFF), 'extraction_parameters/')

    for fiber in old_p_id_DFFFF.keys():
        for order in old_p_id_DFFFF[fiber].keys():
            np.testing.assert_allclose(p_id_DFFFF[fiber][order], old_p_id_DFFFF[fiber][order],
                rtol=1e-5, atol=1e-5)
            


# =====================================================
# HDF5 helper functions
# =====================================================

def load_mat(name, filepath):
    with h5py.File(filepath, mode="r") as h5f:
        return h5f[name][()]

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