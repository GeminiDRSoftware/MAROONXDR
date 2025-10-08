"""
Test suite for bundle export primitives and recipe.
"""
import logging
from copy import deepcopy

import astrodata
import pytest

import maroonx_instruments  # noqa - import is necessary for astrodata
from maroonxdr.maroonx.primitives_maroonx_spectrum import MAROONXSpectrum


@pytest.mark.parametrize('bundle_filename', ['N20241114M3271.fits'])
def test_separate_and_bundle_arms(caplog, bundle_filename):
    """
    Test that separateArmStreams and bundleArmStreams correctly split and
    re-bundle MAROON-X data.

    This test verifies the complete workflow:
    1. Split a bundle into Blue and Red arms (splitBundle)
    2. Separate arms into streams (separateArmStreams)
    3. Re-bundle them back together (bundleArmStreams)
    4. Verify the result matches the original structure

    Parameters
    ----------
    caplog : fixture
        pytest logging fixture
    bundle_filename : str
        Name of the bundle file to test

    Returns
    -------
    None
    """
    caplog.set_level(logging.DEBUG)

    # Load the bundle
    ad_bundle = astrodata.open(bundle_filename)
    original_filename = ad_bundle.filename

    # Split the bundle into separate arm files
    p = MAROONXSpectrum([ad_bundle])
    arm_list = p.splitBundle()

    assert len(arm_list) == 2, 'Should split into 2 arms'

    # Verify split worked correctly
    arms = [ad[0].hdr.get('ARM') for ad in arm_list]
    assert 'BLUE' in arms and 'RED' in arms, 'Should have BLUE and RED arms'

    # All should have same ARCHNAME
    archnames = [ad.phu.get('ARCHNAME') for ad in arm_list]
    assert len(set(archnames)) == 1, 'All arms should have same ARCHNAME'
    assert archnames[0] == original_filename, 'ARCHNAME should be original filename'

    # Separate arms into streams
    p2 = MAROONXSpectrum(arm_list)
    blue_list = p2.separateArmStreams()

    assert 'RED' in p2.streams, 'RED stream should be created'
    red_list = p2.streams['RED']

    assert len(blue_list) == 1, 'Should have 1 blue arm'
    assert len(red_list) == 1, 'Should have 1 red arm'

    # Verify separation worked correctly
    assert 'BLUE' in blue_list[0].tags, 'Blue list should have BLUE tag'
    assert 'RED' in red_list[0].tags, 'Red list should have RED tag'

    # Re-bundle the arms
    p3 = MAROONXSpectrum(blue_list)
    p3.streams['RED'] = red_list
    bundled_list = p3.bundleArmStreams()

    assert len(bundled_list) == 1, 'Should produce one bundle'

    bundle_ad = bundled_list[0]

    # Verify the bundle structure
    assert len(bundle_ad) == 2, 'Bundle should have 2 extensions'
    assert bundle_ad.filename == original_filename, 'Filename should be restored'

    # Verify both arms are present
    bundled_arms = [bundle_ad[i].hdr.get('ARM') for i in bundle_ad.indices]
    assert 'BLUE' in bundled_arms, 'Bundle should contain BLUE arm'
    assert 'RED' in bundled_arms, 'Bundle should contain RED arm'

    # Verify ORIGNAME is set correctly
    assert bundle_ad.phu.get('ORIGNAME') == original_filename, \
        'ORIGNAME should be original filename'

    # Verify ARCHNAME was removed (it's now the filename)
    assert 'ARCHNAME' not in bundle_ad.phu, \
        'ARCHNAME should be removed after bundling'


@pytest.mark.parametrize('bundle_filename', ['N20241114M3271.fits'])
def test_separate_arms_error_handling(caplog, bundle_filename):
    """
    Test error handling in separateArmStreams when ARCHNAME is missing.

    Parameters
    ----------
    caplog : fixture
        pytest logging fixture
    bundle_filename : str
        Name of the bundle file to test

    Returns
    -------
    None
    """
    caplog.set_level(logging.WARNING)

    # Load and split bundle
    ad_bundle = astrodata.open(bundle_filename)
    p = MAROONXSpectrum([ad_bundle])
    arm_list = p.splitBundle()

    # Remove ARCHNAME from one arm
    arm_list[0].phu['ARCHNAME'] = None

    # Try to separate - should log warning
    p2 = MAROONXSpectrum(arm_list)
    p2.separateArmStreams()

    # Check that appropriate warning was logged
    assert any('No BLUE or RED tag found' in r.message or
              'No BLUE or RED tagged files found' in r.message
              for r in caplog.records), \
        'Should warn about missing/unmatched files'


@pytest.mark.parametrize('bundle_filename', ['N20241114M3271.fits'])
def test_bundle_arms_requires_red_stream(caplog, bundle_filename):
    """
    Test that bundleArmStreams raises error when RED stream is missing.

    Parameters
    ----------
    caplog : fixture
        pytest logging fixture
    bundle_filename : str
        Name of the bundle file to test

    Returns
    -------
    None
    """
    # Load and split bundle
    ad_bundle = astrodata.open(bundle_filename)
    p = MAROONXSpectrum([ad_bundle])
    arm_list = p.splitBundle()

    # Get just the blue arm
    blue_arm = [ad for ad in arm_list if 'BLUE' in ad.tags]

    # Try to bundle without running separateArmStreams first
    p2 = MAROONXSpectrum(blue_arm)

    with pytest.raises(ValueError, match='RED stream not found'):
        p2.bundleArmStreams()


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
