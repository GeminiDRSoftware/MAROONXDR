"""Tests for file sorting and stream separation primitives."""

import logging
from copy import deepcopy

import astrodata
import numpy as np
import pytest

from maroonxdr.maroonx.primitives_maroonx_2D import MAROONX


@pytest.mark.parametrize('bundle_filename', ['N20241114M3271.fits'])
def test_splitBundle(caplog, bundle_filename):
    """
    Test that a Bundle is splitted in Blue and Red astrodata objects.

    Parameters
    ----------
    caplog : fixture
    bundle_filename : str

    Returns
    -------
    None
    """
    caplog.set_level(logging.DEBUG)

    ad_bundle = astrodata.open(bundle_filename)

    p = MAROONX([ad_bundle])
    out = p.splitBundle()

    assert len(out) == 2, 'Should return two astrodata objects'

    ad_1, ad_2 = out
    assert ad_1.indices == [0], 'Should have only one indexable header'
    assert ad_2.indices == [0], 'Should have only one indexable header'

    assert ad_1[0].hdr.get('ARM') in ['BLUE', 'RED'], 'ARM should be BLUE or RED'
    assert ad_2[0].hdr.get('ARM') in ['BLUE', 'RED'], 'ARM should be BLUE or RED'
    assert ad_1[0].hdr.get('ARM') != ad_2[0].hdr.get('ARM')

    assert ad_1.phu.get('ORIGNAME') != ad_2.phu.get('ORIGNAME')
    assert ad_1.phu.get('ARCHNAME') == ad_2.phu.get('ARCHNAME')

    assert ad_bundle.tables == ad_1.tables == ad_2.tables


@pytest.mark.parametrize('filename_r', ['20241114T181028Z_DFFFD_r_0002.fits'])
@pytest.mark.parametrize('filename_b', ['20241114T181028Z_DFFFD_b_0008.fits'])
def test_checkArm_collection_and_rejection(caplog, filename_r, filename_b):
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

    ad_red = astrodata.open(filename_r)
    ad_blue = astrodata.open(filename_b)

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


@pytest.mark.parametrize('DFFFD_file', ['20241114T181815Z_DFFFD_b_0008.fits'])
@pytest.mark.parametrize('FDDDF_file', ['20241114T191006Z_DDDDF_b_0007.fits'])
def test_separating_flat_streams(caplog, DFFFD_file, FDDDF_file):
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
    ad_DFFFD = astrodata.open(DFFFD_file)
    ad_FDDDF = astrodata.open(FDDDF_file)
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
def test_combining_flat_streams(caplog, DFFFD_file, FDDDF_file):
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
    ad_DFFFD = astrodata.open(DFFFD_file)
    ad_FDDDF = astrodata.open(FDDDF_file)
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
        adtest[0].data[0],
        np.max([ad_DFFFD.data[0], ad_FDDDF.data[0]], axis=0)
    )
    np.testing.assert_allclose(
        adtest[0].variance[0],
        np.max([p_var[0].variance[0], p_var[1].variance[0]], axis=0)
    )


if __name__ == '__main__':
    pytest.main()
