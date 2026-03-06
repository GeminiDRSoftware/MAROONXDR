"""
Unit tests for splitBundle primitive.

This is a suite of tests to be run with pytest.
"""

import os
import shutil

import astrodata
import pytest
from astrodata.testing import download_from_archive

import maroonx_instruments  # noqa - import is necessary for astrodata
from maroonxdr.maroonx.primitives_maroonx_2D import MAROONX

# -- Test datasets -------------------------------------------------------------
test_datasets = ['N20241114M3271.fits']
exposuremeter_datasets = ['N20241124M1413.fits']


# -- Tests ---------------------------------------------------------------------
@pytest.mark.parametrize('filename', test_datasets)
def test_split_bundle(path_to_inputs, filename):
    """
    Test that splitBundle() correctly splits a MAROON-X bundle into separate arms.
    """
    ad = astrodata.open(os.path.join(path_to_inputs, filename))

    # Verify this is a bundle with 2 extensions
    assert len(ad) == 2, 'Input should be a bundle with 2 extensions'
    assert 'BUNDLE' in ad.tags, 'Input should have BUNDLE tag'

    # Store original filename
    original_filename = ad.filename

    # Run splitBundle
    p = MAROONX([ad])
    split_ads = p.splitBundle()

    # Should produce 2 outputs (one for each arm)
    assert len(split_ads) == 2, 'splitBundle should produce 2 outputs'

    arm_1, arm_2 = split_ads

    # Each output should have only 1 extension
    assert len(arm_1) == 1, 'should have 1 extension'
    assert len(arm_2) == 1, 'should have 1 extension'

    # Collect ARM values
    assert arm_1.arm()[0] == 'BLUE', 'Extension 1 should be BLUE arm'
    assert arm_2.arm()[0] == 'RED', 'Extension 2 should be RED arm'

    # Verify ARCHNAME references the original bundle
    for arm_ad in split_ads:
        assert (
            arm_ad.phu.get('ARCHNAME') == original_filename
        ), 'ARCHNAME should reference original bundle filename'

    # Verify ORIGNAME is set from extension header
    for arm_ad in split_ads:
        origname = arm_ad.phu.get('ORIGNAME')
        assert origname is not None, 'ORIGNAME should be set'
        assert arm_ad.filename == origname, 'Filename should match ORIGNAME'


@pytest.mark.parametrize('filename', exposuremeter_datasets)
def test_splitBundle_exposuremeter(path_to_inputs, filename):
    """
    Test that splitBundle preserves EXPOSUREMETER table metadata.

    Parameters
    ----------
    path_to_inputs : fixture
    filename : str
    """
    filepath = os.path.join(path_to_inputs, filename)
    if not os.path.isfile(filepath):
        pytest.skip(f'{filename} not available for testing.')

    ad_bundle = astrodata.open(filepath)

    p = MAROONX([ad_bundle])
    out = p.splitBundle()

    ad_1, ad_2 = out

    assert (
        ad_1.EXPOSUREMETER.meta['header']['ZP_PC']
        == ad_bundle.EXPOSUREMETER.meta['header']['TZERO2']
    )
    assert (
        ad_1.EXPOSUREMETER.meta['header']['ZP_FRD']
        == ad_bundle.EXPOSUREMETER.meta['header']['TZERO3']
    )


# -- Create inputs -------------------------------------------------------------
def create_inputs():
    """
    Create input files for this test module.

    Run with: python -m maroonxdr.maroonx.tests.bundle.test_bundle --create-inputs
    """
    input_path = os.path.join(
        os.environ['DRAGONS_TEST'],
        'maroonxdr',
        'maroonx',
        'bundle',
        'test_bundle',
        'inputs',
    )
    os.makedirs(input_path, exist_ok=True)

    # Download raw files to default raw_files/ cache
    for filename in test_datasets + exposuremeter_datasets:
        print(f'  Downloading {filename}')
        raw_path = download_from_archive(filename)

        # No preprocessing needed — raw bundles are the test inputs
        shutil.copy2(raw_path, os.path.join(input_path, filename))
        print(f'  Copied to {input_path}')


if __name__ == '__main__':
    import sys

    if '--create-inputs' in sys.argv[1:]:
        create_inputs()
    else:
        pytest.main()
