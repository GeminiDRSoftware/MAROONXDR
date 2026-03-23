import os

import astrodata
import h5py
import numpy as np
import pytest
from astropy.io import fits
from gempy.adlibrary import dataselect

import maroonx_instruments  # noqa : important to load adclass tags
from maroonxdr.maroonx.primitives_maroonx_2D import MAROONX

from . import legacy_adapter

# =========================================================
# TESTS
# =========================================================

@pytest.mark.preprocessed_data
def test_stackFlats(arm, path_to_legacy_flats, preprocessed_files_path):
    """
    Needs legacy output of step 1, create master flats:
    reduce/recipes/make_master_flats.py
    """
    if arm == "BLUE":
        old_DFFFD = path_to_legacy_flats / "20241114T18_masterflat_DFFFD_b_0008.fits"
        old_DDDDF = path_to_legacy_flats / "20241114T19_masterflat_DDDDF_b_0007.fits"
    else:
        old_DFFFD = path_to_legacy_flats / "20241114T18_masterflat_DFFFD_r_0002.fits"
        old_DDDDF = path_to_legacy_flats / "20241114T19_masterflat_DDDDF_r_0002.fits"

    # get all flat files
    raw_files = sorted([str(f) for f in preprocessed_files_path.glob('*.fits')])
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

    # Compare with old masterflats.
    # Legacy stores stacked flats as float64 but casts to float32 downstream
    # (backgroundfit.py line 237). Compare at float32 to match the precision
    # that actually matters for the final product.
    with fits.open(old_DDDDF) as hdu_DDDDF, fits.open(old_DFFFD) as hdu_DFFFD:

        np.testing.assert_allclose(
            ad_DDDDF[0].data.astype(np.float32),
            hdu_DDDDF[0].data.astype(np.float32),
            rtol=0, atol=1e-4)
        np.testing.assert_allclose(
            ad_DFFFD[0].data.astype(np.float32),
            hdu_DFFFD[0].data.astype(np.float32),
            rtol=0, atol=1e-4)


@pytest.mark.preprocessed_data
def test_identifyStripes(arm, path_to_legacy_flats, preprocessed_files_path):
    """
    Needs legacy output of step 1, create master flats:
    reduce/extraction.py
    """
    if arm == "BLUE":
        old_DFFFD = path_to_legacy_flats / "20241114T18_masterflat_DFFFD_b_0008.hdf"
        old_DDDDF = path_to_legacy_flats / "20241114T19_masterflat_DDDDF_b_0007.hdf"
    else:
        old_DFFFD = path_to_legacy_flats / "20241114T18_masterflat_DFFFD_r_0002.hdf"
        old_DDDDF = path_to_legacy_flats / "20241114T19_masterflat_DDDDF_r_0002.hdf"

    # get all flat files
    raw_files = sorted([str(f) for f in preprocessed_files_path.glob('*.fits')])
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

    p.identifyStripes(selected_fibers=[5]) # identify stripes based on MX architecture files
    p.identifyStripes(stream='DFFFD_flats',selected_fibers=[2,3,4])

    p_id_DDDDF = p.streams['main'][0][0].STRIPES_ID
    p_id_DFFFD = p.streams['DFFFD_flats'][0][0].STRIPES_ID

    old_p_id_DDDDF = legacy_adapter.load_dict_from_hdf5(str(old_DDDDF), 'extraction_parameters/')
    old_p_id_DFFFD = legacy_adapter.load_dict_from_hdf5(str(old_DFFFD), 'extraction_parameters/')

    for fiber in old_p_id_DDDDF.keys():
        for order in old_p_id_DDDDF[fiber].keys():
            np.testing.assert_allclose(p_id_DDDDF[fiber][order], old_p_id_DDDDF[fiber][order],
                rtol=0.03, atol=0)

    for fiber in old_p_id_DFFFD.keys():
        for order in old_p_id_DFFFD[fiber].keys():
            np.testing.assert_allclose(p_id_DFFFD[fiber][order], old_p_id_DFFFD[fiber][order],
                rtol=0.03, atol=0)


@pytest.mark.preprocessed_data
def test_identifyStripes_legacyOrder(arm, path_to_legacy_flats, preprocessed_files_path):
    """
    Needs legacy output of step 1, create master flats:
    reduce/extraction.py
    """
    if arm == "BLUE":
        old_DFFFD = path_to_legacy_flats / "20241114T18_masterflat_DFFFD_b_0008.hdf"
        old_DDDDF = path_to_legacy_flats / "20241114T19_masterflat_DDDDF_b_0007.hdf"
    else:
        old_DFFFD = path_to_legacy_flats / "20241114T18_masterflat_DFFFD_r_0002.hdf"
        old_DDDDF = path_to_legacy_flats / "20241114T19_masterflat_DDDDF_r_0002.hdf"

    # get all flat files
    raw_files = sorted([str(f) for f in preprocessed_files_path.glob('*.fits')])
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

    p.identifyStripes(stream='main', selected_fibers=[5])   # '0,0,0,0,5') # identify stripes based on MX architecture files
    p.identifyStripes(stream='DFFFD_flats', selected_fibers=[2, 3, 4]) # '0,2,3,4,0')

    p_id_DDDDF = p.streams['main'][0][0].STRIPES_ID
    p_id_DFFFD = p.streams['DFFFD_flats'][0][0].STRIPES_ID

    old_p_id_DDDDF = legacy_adapter.load_dict_from_hdf5(str(old_DDDDF), 'extraction_parameters/')
    old_p_id_DFFFD = legacy_adapter.load_dict_from_hdf5(str(old_DFFFD), 'extraction_parameters/')

    for fiber in old_p_id_DDDDF.keys():
        for order in old_p_id_DDDDF[fiber].keys():
            np.testing.assert_allclose(p_id_DDDDF[fiber][order], old_p_id_DDDDF[fiber][order],
                rtol=0.03, atol=0)

    for fiber in old_p_id_DFFFD.keys():
        for order in old_p_id_DFFFD[fiber].keys():
            np.testing.assert_allclose(p_id_DFFFD[fiber][order], old_p_id_DFFFD[fiber][order],
                rtol=0.03, atol=0)


@pytest.mark.preprocessed_data
@pytest.mark.xfail(
    reason="Uses single overscan flow; legacy .npy snapshots were recorded "
           "with double overscan (backgroundfit.py). Use "
           "test_removeStraylight_legacyOrder for the correct legacy flow.",
    strict=True,
)
def test_removeStraylight(arm, path_to_legacy_flats, preprocessed_files_path):
    # This test does not subtract overscan twice as test_removeStraylight_legacyOrder does.

    if arm == "BLUE":
        old_DFFFD = path_to_legacy_flats / "20241114T18_masterflat_backgroundsubtracted_DFFFD_b_0008.fits"
        old_DDDDF = path_to_legacy_flats / "20241114T19_masterflat_backgroundsubtracted_DDDDF_b_0007.fits"
    else:
        old_DFFFD = path_to_legacy_flats / "20241114T18_masterflat_backgroundsubtracted_DFFFD_r_0002.fits"
        old_DDDDF = path_to_legacy_flats / "20241114T19_masterflat_backgroundsubtracted_DDDDF_r_0002.fits"

    # get all flat files
    raw_files = sorted([str(f) for f in preprocessed_files_path.glob('*.fits')])
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

    p.identifyStripes(selected_fibers=[5]) # identify stripes based on MX architecture files
    p.identifyStripes(stream='DFFFD_flats',selected_fibers=[2,3,4])

    p.defineFlatStripes()  # defines pixel inclusion for each flat region based on stripe ids
    p.defineFlatStripes(stream='DFFFD_flats')

    p.removeStrayLight(filter_size=19, box_size=20, snapshot=False, legacy=True)  # remove straylight from frame (this is why 2 partial illumination flat sets are necessary)
    p.removeStrayLight(stream='DFFFD_flats', filter_size=19, box_size=20, snapshot=False, legacy=True)

    ad_DDDDF = p.streams['main'][0]
    ad_DFFFD = p.streams['DFFFD_flats'][0]

    # Compare with old masterflats
    with fits.open(old_DDDDF) as hdu_DDDDF, fits.open(old_DFFFD) as hdu_DFFFD:

        np.testing.assert_allclose(
            ad_DDDDF[0].data.astype(np.float32),
            hdu_DDDDF[0].data.astype(np.float32),
            rtol=0, atol=1e-4)
        np.testing.assert_allclose(
            ad_DFFFD[0].data.astype(np.float32),
            hdu_DFFFD[0].data.astype(np.float32),
            rtol=0, atol=1e-4)


@pytest.mark.preprocessed_data
def test_removeStraylight_legacyOrder(arm, path_to_legacy_flats, preprocessed_files_path):

    if arm == "BLUE":
        old_DFFFD = path_to_legacy_flats / "20241114T18_masterflat_backgroundsubtracted_DFFFD_b_0008.fits"
        old_DDDDF = path_to_legacy_flats / "20241114T19_masterflat_backgroundsubtracted_DDDDF_b_0007.fits"
    else:
        old_DFFFD = path_to_legacy_flats / "20241114T18_masterflat_backgroundsubtracted_DFFFD_r_0002.fits"
        old_DDDDF = path_to_legacy_flats / "20241114T19_masterflat_backgroundsubtracted_DDDDF_r_0002.fits"

    # get all flat files
    raw_files = sorted([str(f) for f in preprocessed_files_path.glob('*.fits')])
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

    p.identifyStripes(selected_fibers=[5]) # identify stripes based on MX architecture files
    p.identifyStripes(stream='DFFFD_flats',selected_fibers=[2,3,4])

    p.defineFlatStripes()  # defines pixel inclusion for each flat region based on stripe ids
    p.defineFlatStripes(stream='DFFFD_flats')

    p.removeStrayLight(filter_size=19, box_size=20, snapshot=False, legacy=True)  # remove straylight from frame (this is why 2 partial illumination flat sets are necessary)
    p.removeStrayLight(stream='DFFFD_flats', filter_size=19, box_size=20, snapshot=False, legacy=True)

    ad_DDDDF = p.streams['main'][0]
    ad_DFFFD = p.streams['DFFFD_flats'][0]

    # Compare with old masterflats.
    # Legacy backgroundfit.py processes at float32 (line 237) but saves as
    # float64. DRAGONS legacy patch also produces float32. Compare at float32.
    with fits.open(old_DDDDF) as hdu_DDDDF, fits.open(old_DFFFD) as hdu_DFFFD:

        np.testing.assert_allclose(
            ad_DDDDF[0].data.astype(np.float32),
            hdu_DDDDF[0].data.astype(np.float32),
            rtol=0, atol=1e-4)
        np.testing.assert_allclose(
            ad_DFFFD[0].data.astype(np.float32),
            hdu_DFFFD[0].data.astype(np.float32),
            rtol=0, atol=1e-4)


@pytest.mark.preprocessed_data
def test_combinedFlat(arm, path_to_legacy_flats, preprocessed_files_path):

    if arm == "BLUE":
        old_DFFFF = path_to_legacy_flats / "20241114T19_masterflat_backgroundsubtracted_FFFFF_b_0007.hdf"
    else:
        old_DFFFF = path_to_legacy_flats / "20241114T19_masterflat_backgroundsubtracted_FFFFF_r_0002.hdf"

    # get all flat files
    raw_files = sorted([str(f) for f in preprocessed_files_path.glob('*.fits')])
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

    p.identifyStripes(selected_fibers=[5]) # identify stripes based on MX architecture files
    p.identifyStripes(stream='DFFFD_flats',selected_fibers=[2,3,4])

    p.defineFlatStripes()  # defines pixel inclusion for each flat region based on stripe ids
    p.defineFlatStripes(stream='DFFFD_flats')

    p.removeStrayLight(filter_size=19, box_size=20, snapshot=False, legacy=True)  # remove straylight from frame (this is why 2 partial illumination flat sets are necessary)
    p.removeStrayLight(stream='DFFFD_flats', filter_size=19, box_size=20, snapshot=False, legacy=True)
    # =========================================================

    p.combineFlatStreams(stream='main', stream_2='DFFFD_flats')  # combine straylight-removed images
    p.clearStream(stream='DFFFD_flats') # remove second stream

    p.findStripes()  # re-run find/identify/define routine on combined frame
    p.identifyStripes(selected_fibers=[2,3,4,5])
    # p.defineFlatStripes(extract=True)

    p_id_DFFFF = p.streams['main'][0][0].STRIPES_ID

    old_p_id_DFFFF = legacy_adapter.load_dict_from_hdf5(str(old_DFFFF), 'extraction_parameters/')

    for fiber in old_p_id_DFFFF.keys():
        for order in old_p_id_DFFFF[fiber].keys():
            np.testing.assert_allclose(p_id_DFFFF[fiber][order], old_p_id_DFFFF[fiber][order],
                rtol=1e-5, atol=1e-5)


# =========================================================
# NEW TESTS
# =========================================================


@pytest.mark.slow
@pytest.mark.preprocessed_data
@pytest.mark.parametrize(
    'matching_filenames',
    [
        (
            '20241114T190714Z_DDDDF_b_0007_DFFFF_flat.fits',
            '20241114T19_masterflat_backgroundsubtracted_FFFFF_b_0007.hdf',
        ),
        (
            '20241114T190714Z_DDDDF_r_0002_DFFFF_flat.fits',
            '20241114T19_masterflat_backgroundsubtracted_FFFFF_r_0002.hdf',
        ),
    ],
)
def test_masterflat_index_maps(path_to_inputs, path_to_legacy_flats, matching_filenames):
    """Compare INDEX_FIBER and INDEX_ORDER pixel maps against legacy stripe_indices."""
    file_name, legacy_file_name = matching_filenames

    ad = astrodata.open(os.path.join(path_to_inputs, file_name))
    legacy_file = str(path_to_legacy_flats / legacy_file_name)

    new_index_fiber = ad[0].INDEX_FIBER
    new_index_order = ad[0].INDEX_ORDER

    old_index_fiber = legacy_adapter.load_mat('stripe_indices/fiber', legacy_file)
    old_index_order = legacy_adapter.load_mat('stripe_indices/order', legacy_file)

    assert new_index_fiber.shape == old_index_fiber.shape
    assert new_index_order.shape == old_index_order.shape

    np.testing.assert_array_equal(new_index_fiber, old_index_fiber,
        err_msg='INDEX_FIBER mismatch')
    np.testing.assert_array_equal(new_index_order, old_index_order,
        err_msg='INDEX_ORDER mismatch')


@pytest.mark.slow
@pytest.mark.preprocessed_data
@pytest.mark.parametrize(
    'matching_filenames',
    [
        (
            '20241114T190714Z_DDDDF_b_0007_DFFFF_flat.fits',
            '20241114T19_masterflat_backgroundsubtracted_FFFFF_b_0007.hdf',
        ),
        (
            '20241114T190714Z_DDDDF_r_0002_DFFFF_flat.fits',
            '20241114T19_masterflat_backgroundsubtracted_FFFFF_r_0002.hdf',
        ),
    ],
)
def test_masterflat_stripes_id(path_to_inputs, path_to_legacy_flats, matching_filenames):
    """Compare STRIPES_ID polynomial coefficients against legacy extraction_parameters."""
    file_name, legacy_file_name = matching_filenames

    ad = astrodata.open(os.path.join(path_to_inputs, file_name))
    legacy_file = str(path_to_legacy_flats / legacy_file_name)

    stripes_table = ad[0].STRIPES_ID
    stripes_fibers = ad[0].STRIPES_FIBERS

    old_p_id = legacy_adapter.load_dict_from_hdf5(legacy_file, 'extraction_parameters/')

    # STRIPES_ID table is a vstack of per-fiber tables.
    # Each fiber contributes n_coeffs rows (polynomial degree + 1).
    n_fibers = len(stripes_fibers)
    n_coeffs = len(stripes_table) // n_fibers

    for i, fiber_num in enumerate(stripes_fibers):
        fiber_key = f'fiber_{fiber_num}'
        start = i * n_coeffs
        end = start + n_coeffs

        for order in old_p_id[fiber_key].keys():
            np.testing.assert_allclose(
                stripes_table[order][start:end], old_p_id[fiber_key][order],
                rtol=1e-5, atol=1e-5,
            )