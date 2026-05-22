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

# Legacy glob() ordering extracted from masterflat FITS HISTORY headers.
# DDDDF files first, then DFFFD — separateFlatStreams preserves within-group order.
LEGACY_FILE_ORDER = {
    'BLUE': [
        '20250701T172509Z_DDDDF_b_0007.fits',
        '20250701T172324Z_DDDDF_b_0007.fits',
        '20250701T172140Z_DDDDF_b_0007.fits',
        '20250701T171955Z_DDDDF_b_0007.fits',
        '20250701T171811Z_DDDDF_b_0007.fits',
        '20250701T171553Z_DDDDF_b_0007.fits',
        '20250701T170537Z_DFFFD_b_0008.fits',
        '20250701T170101Z_DFFFD_b_0008.fits',
        '20250701T171051Z_DFFFD_b_0008.fits',
        '20250701T170906Z_DFFFD_b_0008.fits',
        '20250701T170353Z_DFFFD_b_0008.fits',
        '20250701T170721Z_DFFFD_b_0008.fits',
    ],
    'RED': [
        '20250701T171955Z_DDDDF_r_0002.fits',
        '20250701T172140Z_DDDDF_r_0002.fits',
        '20250701T171553Z_DDDDF_r_0002.fits',
        '20250701T172324Z_DDDDF_r_0002.fits',
        '20250701T171811Z_DDDDF_r_0002.fits',
        '20250701T172509Z_DDDDF_r_0002.fits',
        '20250701T170906Z_DFFFD_r_0002.fits',
        '20250701T170721Z_DFFFD_r_0002.fits',
        '20250701T170101Z_DFFFD_r_0002.fits',
        '20250701T170353Z_DFFFD_r_0002.fits',
        '20250701T171051Z_DFFFD_r_0002.fits',
        '20250701T170537Z_DFFFD_r_0002.fits',
    ],
}


def _legacy_ordered_files(preprocessed_files_path, arm):
    """Return full paths to flat files in the legacy glob() order."""
    return [str(preprocessed_files_path / f) for f in LEGACY_FILE_ORDER[arm]]


@pytest.mark.preprocessed_data
def test_stackFlats(arm, path_to_legacy_flats, preprocessed_files_path):

    if arm == "BLUE":
        old_DFFFD = path_to_legacy_flats / "20250701T17_masterflat_DFFFD_b_0008.fits"
        old_DDDDF = path_to_legacy_flats / "20250701T17_masterflat_DDDDF_b_0007.fits"
    else:
        old_DFFFD = path_to_legacy_flats / "20250701T17_masterflat_DFFFD_r_0002.fits"
        old_DDDDF = path_to_legacy_flats / "20250701T17_masterflat_DDDDF_r_0002.fits"

    selected_spect = _legacy_ordered_files(preprocessed_files_path, arm)

    adinput = [astrodata.open(f) for f in selected_spect]
    p = MAROONX(adinput)

    p.subtractOverscan()
    p.separateFlatStreams()

    p.stackFlats(stream='main', scale_mode='first_frame', suffix='_DDDDF_flat')
    p.stackFlats(stream='DFFFD_flats', scale_mode='first_frame', suffix='_DFFFD_flat')

    ad_DDDDF = p.streams['main'][0]
    ad_DFFFD = p.streams['DFFFD_flats'][0]

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

    if arm == "BLUE":
        old_DFFFD = path_to_legacy_flats / "20250701T17_masterflat_DFFFD_b_0008.hdf"
        old_DDDDF = path_to_legacy_flats / "20250701T17_masterflat_DDDDF_b_0007.hdf"
    else:
        old_DFFFD = path_to_legacy_flats / "20250701T17_masterflat_DFFFD_r_0002.hdf"
        old_DDDDF = path_to_legacy_flats / "20250701T17_masterflat_DDDDF_r_0002.hdf"

    selected_spect = _legacy_ordered_files(preprocessed_files_path, arm)

    adinput = [astrodata.open(f) for f in selected_spect]
    p = MAROONX(adinput)

    p.subtractOverscan()
    p.trimOverscan()
    p.correctImageOrientation()
    p.separateFlatStreams()

    p.stackFlats(stream='main', scale_mode='first_frame', suffix='_DDDDF_flat')
    p.stackFlats(stream='DFFFD_flats', scale_mode='first_frame', suffix='_DFFFD_flat')

    p.findStripes()
    p.findStripes(stream='DFFFD_flats')

    p.identifyStripes(selected_fibers=[5])
    p.identifyStripes(stream='DFFFD_flats', selected_fibers=[2, 3, 4])

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

    if arm == "BLUE":
        old_DFFFD = path_to_legacy_flats / "20250701T17_masterflat_DFFFD_b_0008.hdf"
        old_DDDDF = path_to_legacy_flats / "20250701T17_masterflat_DDDDF_b_0007.hdf"
    else:
        old_DFFFD = path_to_legacy_flats / "20250701T17_masterflat_DFFFD_r_0002.hdf"
        old_DDDDF = path_to_legacy_flats / "20250701T17_masterflat_DDDDF_r_0002.hdf"

    selected_spect = _legacy_ordered_files(preprocessed_files_path, arm)

    adinput = [astrodata.open(f) for f in selected_spect]
    p = MAROONX(adinput)

    p.prepare()
    p.checkArm()
    p.checkND()
    p.addDQ()
    p.addVAR(read_noise=True, poisson_noise=True)

    p.subtractOverscan()

    p.separateFlatStreams()

    p.stackFlats(stream='main', scale_mode='first_frame', suffix='_DDDDF_flat')
    p.stackFlats(stream='DFFFD_flats', scale_mode='first_frame', suffix='_DFFFD_flat')

    p.subtractOverscan(stream='main')
    p.subtractOverscan(stream='DFFFD_flats')

    p.trimOverscan(stream='main')
    p.trimOverscan(stream='DFFFD_flats')

    p.correctImageOrientation(stream='main')
    p.correctImageOrientation(stream='DFFFD_flats')

    p.findStripes(stream='main')
    p.findStripes(stream='DFFFD_flats')

    p.identifyStripes(stream='main', selected_fibers=[5])
    p.identifyStripes(stream='DFFFD_flats', selected_fibers=[2, 3, 4])

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

    if arm == "BLUE":
        old_DFFFD = path_to_legacy_flats / "20250701T17_masterflat_backgroundsubtracted_DFFFD_b_0008.fits"
        old_DDDDF = path_to_legacy_flats / "20250701T17_masterflat_backgroundsubtracted_DDDDF_b_0007.fits"
    else:
        old_DFFFD = path_to_legacy_flats / "20250701T17_masterflat_backgroundsubtracted_DFFFD_r_0002.fits"
        old_DDDDF = path_to_legacy_flats / "20250701T17_masterflat_backgroundsubtracted_DDDDF_r_0002.fits"

    selected_spect = _legacy_ordered_files(preprocessed_files_path, arm)

    adinput = [astrodata.open(f) for f in selected_spect]
    p = MAROONX(adinput)

    p.prepare()
    p.checkArm()
    p.checkND()
    p.addDQ()
    p.addVAR(read_noise=True, poisson_noise=True)
    p.subtractOverscan()
    p.trimOverscan()
    p.correctImageOrientation()
    p.separateFlatStreams()

    p.stackFlats(stream='main', scale_mode='first_frame', suffix='_flat')
    p.stackFlats(stream='DFFFD_flats', scale_mode='first_frame', suffix='_flat')

    p.findStripes()
    p.findStripes(stream='DFFFD_flats')

    p.identifyStripes(selected_fibers=[5])
    p.identifyStripes(stream='DFFFD_flats', selected_fibers=[2, 3, 4])

    p.defineFlatStripes()
    p.defineFlatStripes(stream='DFFFD_flats')

    p.removeStrayLight(filter_size=19, box_size=20, snapshot=False, legacy=True)
    p.removeStrayLight(stream='DFFFD_flats', filter_size=19, box_size=20, snapshot=False, legacy=True)

    ad_DDDDF = p.streams['main'][0]
    ad_DFFFD = p.streams['DFFFD_flats'][0]

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
        old_DFFFD = path_to_legacy_flats / "20250701T17_masterflat_backgroundsubtracted_DFFFD_b_0008.fits"
        old_DDDDF = path_to_legacy_flats / "20250701T17_masterflat_backgroundsubtracted_DDDDF_b_0007.fits"
    else:
        old_DFFFD = path_to_legacy_flats / "20250701T17_masterflat_backgroundsubtracted_DFFFD_r_0002.fits"
        old_DDDDF = path_to_legacy_flats / "20250701T17_masterflat_backgroundsubtracted_DDDDF_r_0002.fits"

    selected_spect = _legacy_ordered_files(preprocessed_files_path, arm)

    adinput = [astrodata.open(f) for f in selected_spect]
    p = MAROONX(adinput)

    p.prepare()
    p.checkArm()
    p.checkND()
    p.addDQ()
    p.addVAR(read_noise=True, poisson_noise=True)
    p.subtractOverscan()

    p.separateFlatStreams()

    p.stackFlats(stream='main', scale_mode='first_frame', suffix='_flat')
    p.stackFlats(stream='DFFFD_flats', scale_mode='first_frame', suffix='_flat')

    p.subtractOverscan(stream='main')
    p.subtractOverscan(stream='DFFFD_flats')

    p.trimOverscan(stream='main')
    p.trimOverscan(stream='DFFFD_flats')

    p.correctImageOrientation(stream='main')
    p.correctImageOrientation(stream='DFFFD_flats')

    p.findStripes()
    p.findStripes(stream='DFFFD_flats')

    p.identifyStripes(selected_fibers=[5])
    p.identifyStripes(stream='DFFFD_flats', selected_fibers=[2, 3, 4])

    p.defineFlatStripes()
    p.defineFlatStripes(stream='DFFFD_flats')

    p.removeStrayLight(filter_size=19, box_size=20, snapshot=False, legacy=True)
    p.removeStrayLight(stream='DFFFD_flats', filter_size=19, box_size=20, snapshot=False, legacy=True)

    ad_DDDDF = p.streams['main'][0]
    ad_DFFFD = p.streams['DFFFD_flats'][0]

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
        old_FFFFF = path_to_legacy_flats / "20250701T17_masterflat_backgroundsubtracted_FFFFF_b_0007.hdf"
    else:
        old_FFFFF = path_to_legacy_flats / "20250701T17_masterflat_backgroundsubtracted_FFFFF_r_0002.hdf"

    selected_spect = _legacy_ordered_files(preprocessed_files_path, arm)

    adinput = [astrodata.open(f) for f in selected_spect]
    p = MAROONX(adinput)

    p.prepare()
    p.checkArm()
    p.checkND()
    p.addDQ()
    p.addVAR(read_noise=True, poisson_noise=True)
    p.subtractOverscan()
    p.separateFlatStreams()

    p.stackFlats(stream='main', scale_mode='first_frame', suffix='_flat')
    p.stackFlats(stream='DFFFD_flats', scale_mode='first_frame', suffix='_flat')

    p.subtractOverscan(stream='main')
    p.subtractOverscan(stream='DFFFD_flats')

    p.trimOverscan(stream='main')
    p.trimOverscan(stream='DFFFD_flats')

    p.correctImageOrientation(stream='main')
    p.correctImageOrientation(stream='DFFFD_flats')

    p.findStripes()
    p.findStripes(stream='DFFFD_flats')

    p.identifyStripes(selected_fibers=[5])
    p.identifyStripes(stream='DFFFD_flats', selected_fibers=[2, 3, 4])

    p.defineFlatStripes()
    p.defineFlatStripes(stream='DFFFD_flats')

    p.removeStrayLight(filter_size=19, box_size=20, snapshot=False, legacy=True)
    p.removeStrayLight(stream='DFFFD_flats', filter_size=19, box_size=20, snapshot=False, legacy=True)

    p.combineFlatStreams(stream='main', stream_2='DFFFD_flats')
    p.clearStream(stream='DFFFD_flats')

    p.findStripes()
    p.identifyStripes(selected_fibers=[2, 3, 4, 5])

    p_id_FFFFF = p.streams['main'][0][0].STRIPES_ID

    old_p_id_FFFFF = legacy_adapter.load_dict_from_hdf5(str(old_FFFFF), 'extraction_parameters/')

    for fiber in old_p_id_FFFFF.keys():
        for order in old_p_id_FFFFF[fiber].keys():
            np.testing.assert_allclose(p_id_FFFFF[fiber][order], old_p_id_FFFFF[fiber][order],
                rtol=1e-5, atol=1e-5)


@pytest.mark.slow
@pytest.mark.preprocessed_data
@pytest.mark.parametrize(
    'matching_filenames',
    [
        (
            '20250701T172509Z_DDDDF_b_0007_DFFFF_flat.fits',
            '20250701T17_masterflat_backgroundsubtracted_FFFFF_b_0007.hdf',
        ),
        (
            '20250701T171955Z_DDDDF_r_0002_DFFFF_flat.fits',
            '20250701T17_masterflat_backgroundsubtracted_FFFFF_r_0002.hdf',
        ),
    ],
)
def test_masterflat_index_maps(path_to_inputs, path_to_legacy_flats, matching_filenames):
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
            '20250701T172509Z_DDDDF_b_0007_DFFFF_flat.fits',
            '20250701T17_masterflat_backgroundsubtracted_FFFFF_b_0007.hdf',
        ),
        (
            '20250701T171955Z_DDDDF_r_0002_DFFFF_flat.fits',
            '20250701T17_masterflat_backgroundsubtracted_FFFFF_r_0002.hdf',
        ),
    ],
)
def test_masterflat_stripes_id(path_to_inputs, path_to_legacy_flats, matching_filenames):
    file_name, legacy_file_name = matching_filenames

    ad = astrodata.open(os.path.join(path_to_inputs, file_name))
    legacy_file = str(path_to_legacy_flats / legacy_file_name)

    stripes_table = ad[0].STRIPES_ID
    stripes_fibers = ad[0].STRIPES_FIBERS

    old_p_id = legacy_adapter.load_dict_from_hdf5(legacy_file, 'extraction_parameters/')

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
