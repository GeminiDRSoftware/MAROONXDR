"""Tests for neutral density filter checking primitives."""

import logging
import os
from copy import deepcopy

import astrodata
import pytest

import maroonx_instruments  # noqa - import is necessary for astrodata
from maroonxdr.maroonx.primitives_maroonx_2D import MAROONX

# -- Test datasets -------------------------------------------------------------
# These bundles are needed for debundling into split-arm files
bundles_needed = [
    'N20241114M3271.fits',
]


# -- Tests ---------------------------------------------------------------------
@pytest.mark.parametrize('filename', ['20241114T181028Z_DFFFD_r_0002.fits'])
def test_nd_filter_good_series(caplog, path_to_inputs, filename):
    """Test ND filter checker outputs full set when all inputs agree with first file.

    Parameters
    ----------
    caplog : fixture
    filename : str

    Returns
    -------
    None
    """
    caplog.set_level(logging.DEBUG)

    ad = astrodata.open(os.path.join(path_to_inputs, filename))
    test_objects = [deepcopy(ad), deepcopy(ad), deepcopy(ad)]
    p = MAROONX(test_objects)
    p.prepare()
    adtest = p.checkND()

    assert len(caplog.records) > 0
    assert len(adtest) == len(test_objects)


@pytest.mark.parametrize('filename', ['20241114T181028Z_DFFFD_r_0002.fits'])
def test_nd_filter_subgood_series(caplog, path_to_inputs, filename):
    """Test ND filter checker outputs subset matching first file's ND filter value.

    Parameters
    ----------
    caplog : fixture
    filename : str

    Returns
    -------
    None
    """
    caplog.set_level(logging.DEBUG)

    ad = astrodata.open(os.path.join(path_to_inputs, filename))
    ad_1 = deepcopy(ad)
    ad_1[0].hdr['HIERARCH MAROONX ND POSITION'] -= 12.3

    test_objects = [deepcopy(ad), deepcopy(ad), deepcopy(ad_1), deepcopy(ad_1)]
    p = MAROONX(test_objects)
    p.prepare()
    adtest = p.checkND()

    assert len(caplog.records) > 0
    assert any('Not all frames have ' in r.message for r in caplog.records)
    assert len(adtest) == sum(
        [
            test_objects[0].filter_orientation()['ND'] == ad.filter_orientation()['ND']
            for ad in test_objects
        ]
    )
    assert all(
        test_objects[-1].filter_orientation()['ND'] != ad.filter_orientation()['ND']
        for ad in adtest
    )


@pytest.mark.parametrize('filename', ['20241114T181028Z_DFFFD_r_0002.fits'])
def test_nd_filter_bad_series(caplog, path_to_inputs, filename):
    """Test ND filter checker raises OSError when only first file has unique ND setting.

    Parameters
    ----------
    caplog : fixture
    filename : str
    """
    caplog.set_level(logging.DEBUG)

    ad = astrodata.open(os.path.join(path_to_inputs, filename))
    ad_1 = deepcopy(ad)
    ad_1[0].hdr['HIERARCH MAROONX ND POSITION'] -= 12.3

    test_objects = [deepcopy(ad), deepcopy(ad_1), deepcopy(ad_1), deepcopy(ad_1)]
    p = MAROONX(test_objects)
    p.prepare()

    with pytest.raises(OSError):
        p.checkND()

    assert any(
        'Only first frame found, of given, with its simcal ND filter setting'
        in r.message
        for r in caplog.records
    )


# -- Create inputs -------------------------------------------------------------
def create_inputs():
    """
    Create input files for this test module.

    Run with: python -m maroonxdr.maroonx.tests.image.test_ND_filter_check --create-inputs
    """
    from astrodata.testing import download_from_archive

    input_path = os.path.join(
        os.environ['DRAGONS_TEST'],
        'maroonxdr',
        'maroonx',
        'image',
        'test_ND_filter_check',
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
