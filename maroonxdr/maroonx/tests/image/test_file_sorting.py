"""Tests for file sorting and stream separation primitives."""

import logging
import os
from copy import deepcopy

import astrodata
import numpy as np
import pytest

import maroonx_instruments  # noqa - import is necessary for astrodata
from maroonxdr.maroonx.primitives_maroonx_2D import MAROONX

# -- Test datasets -------------------------------------------------------------
# These bundles are needed for debundling into split-arm files
bundles_needed = [
    'N20241114M3271.fits',
    'N20241114M3295.fits',
    'N20241114M3450.fits',
]


# -- Tests ---------------------------------------------------------------------
@pytest.mark.parametrize('filename_r', ['20241114T181028Z_DFFFD_r_0002.fits'])
@pytest.mark.parametrize('filename_b', ['20241114T181028Z_DFFFD_b_0008.fits'])
def test_checkArm_collection_and_rejection(
    caplog, path_to_inputs, filename_r, filename_b
):
    """Test that first file and others of its arm-type are included while else are removed.

    Parameters
    ----------
    caplog : fixture
    filename_r : str
    filename_b : str

    Returns
    -------
    None
    """
    caplog.set_level(logging.DEBUG)

    ad_red = astrodata.open(os.path.join(path_to_inputs, filename_r))
    ad_blue = astrodata.open(os.path.join(path_to_inputs, filename_b))

    test_objects = [
        ad_red,
        deepcopy(ad_red),
        deepcopy(ad_red),
        ad_blue,
        deepcopy(ad_blue),
    ]

    p = MAROONX(test_objects)
    p.prepare()
    out = p.checkArm()

    assert len(caplog.records) > 0
    assert any(
        'Not all frames taken with the same camera arm' in r.message
        for r in caplog.records
    )
    assert all(test_objects[0].filename in ad.filename for ad in out)
    assert len(out) == sum(
        [test_objects[0].filename in ad.filename for ad in test_objects]
    )
    assert all(test_objects[-1].filename not in ad.filename for ad in out)


def test_checkArm(ad_min):
    """Mismatched arm frame is dropped, same arm frames are kept."""

    ad1 = ad_min
    ad2 = deepcopy(ad1)
    ad3 = deepcopy(ad1)
    other_arm = 'BLUE' if 'RED' in ad1.tags else 'RED'
    ad3[0].hdr['ARM'] = other_arm

    p = MAROONX([])
    out = p.checkArm([ad1, ad2, ad3])

    assert len(out) == 2
    assert all(ad.filename == ad1.filename for ad in out)


@pytest.mark.parametrize('DFFFD_file', ['20241114T181815Z_DFFFD_b_0008.fits'])
@pytest.mark.parametrize('FDDDF_file', ['20241114T191006Z_DDDDF_b_0007.fits'])
def test_separating_flat_streams(caplog, path_to_inputs, DFFFD_file, FDDDF_file):
    """Test that seperateFlatStreams correctly separates flats by illuminated fibers.

    Parameters
    ----------
    caplog : fixture
    DFFFD_file : str
    FDDDF_file : str

    Returns
    -------
    None
    """
    caplog.set_level(logging.DEBUG)
    ad_DFFFD = astrodata.open(os.path.join(path_to_inputs, DFFFD_file))
    ad_FDDDF = astrodata.open(os.path.join(path_to_inputs, FDDDF_file))
    test_flats = [
        deepcopy(ad_FDDDF),
        deepcopy(ad_FDDDF),
        deepcopy(ad_FDDDF),
        deepcopy(ad_DFFFD),
        deepcopy(ad_DFFFD),
    ]
    p = MAROONX(test_flats)
    p.prepare()
    p.separateFlatStreams()

    assert not any('Not registered as Flat' in r.message for r in caplog.records)
    assert not any('No FDDDF Flats in input list' in r.message for r in caplog.records)
    assert not any('No DFFFD Flats in input list' in r.message for r in caplog.records)

    DARK, FLAT = 'Dark', 'Flat lamp'
    num_DDDDF = sum(
        [ad.fiber_setup() == [DARK, DARK, DARK, DARK, FLAT] for ad in test_flats]
    )
    num_FDDDF = sum(
        [ad.fiber_setup() == [FLAT, DARK, DARK, DARK, FLAT] for ad in test_flats]
    )
    num_DFFFD = sum(
        [ad.fiber_setup() == [DARK, FLAT, FLAT, FLAT, DARK] for ad in test_flats]
    )

    assert len(p.streams['main']) == num_DDDDF + num_FDDDF

    assert p.streams['DFFFD_flats'] is not None
    assert len(p.streams['DFFFD_flats']) == num_DFFFD


@pytest.mark.parametrize('DFFFD_file', ['20241114T181815Z_DFFFD_b_0008.fits'])
@pytest.mark.parametrize('FDDDF_file', ['20241114T191006Z_DDDDF_b_0007.fits'])
def test_combining_flat_streams(caplog, path_to_inputs, DFFFD_file, FDDDF_file):
    """Test that combineFlatStreams creates by-pixel max of DFFFD and FDDDF flats.

    Parameters
    ----------
    caplog : fixture
    DFFFD_file : str
    FDDDF_file : str

    Returns
    -------
    None
    """
    caplog.set_level(logging.DEBUG)
    ad_DFFFD = astrodata.open(os.path.join(path_to_inputs, DFFFD_file))
    ad_FDDDF = astrodata.open(os.path.join(path_to_inputs, FDDDF_file))
    test_flats = [deepcopy(ad_FDDDF), deepcopy(ad_DFFFD)]
    p = MAROONX(test_flats)
    p.prepare()
    p_var = p.addVAR()
    p.separateFlatStreams()
    adtest = p.combineFlatStreams(stream_2='DFFFD_flats')

    assert not any(
        'does not exist so nothing to transfer' in r.message for r in caplog.records
    )
    assert not any('Unexpected stream lengths' in r.message for r in caplog.records)
    assert len(adtest) == 1

    np.testing.assert_array_equal(
        adtest[0].data[0], np.max([ad_DFFFD.data[0], ad_FDDDF.data[0]], axis=0)
    )
    np.testing.assert_allclose(
        adtest[0].variance[0],
        np.max([p_var[0].variance[0], p_var[1].variance[0]], axis=0),
    )


# -- Create inputs -------------------------------------------------------------
def create_inputs():
    """
    Create input files for this test module.

    Run with: python -m maroonxdr.maroonx.tests.image.test_file_sorting --create-inputs
    """
    from astrodata.testing import download_from_archive

    input_path = os.path.join(
        os.environ['DRAGONS_TEST'],
        'maroonxdr',
        'maroonx',
        'image',
        'test_file_sorting',
        'inputs',
    )
    os.makedirs(input_path, exist_ok=True)

    # Download bundles and run splitBundle to produce debundled files
    for filename in bundles_needed:
        print(f'  Downloading {filename}')
        raw_path = download_from_archive(filename)
        ad = astrodata.open(raw_path)
        p = MAROONX([ad])
        split_ads = p.splitBundle()
        for ad_arm in split_ads:
            ad_arm.write(os.path.join(input_path, ad_arm.filename), overwrite=True)
            print(f'  Wrote {ad_arm.filename} to {input_path}')


if __name__ == '__main__':
    import sys

    if '--create-inputs' in sys.argv[1:]:
        create_inputs()
    else:
        pytest.main()
