"""Unit tests for the bundleArmStreams primitive."""

import pytest

from maroonxdr.maroonx.primitives_maroonx_spectrum import MaroonXSpectrum

from . import make_arm

# -- Fixtures ------------------------------------------------------------------
ARCHNAME = 'N00000000M0000.fits'


@pytest.fixture()
def ad_blue():
    """A BLUE arm."""
    return make_arm('BLUE', ARCHNAME)


@pytest.fixture()
def ad_red():
    """A RED arm sharing the ARCHNAME of the blue one."""
    return make_arm('RED', ARCHNAME)


# -- Tests ---------------------------------------------------------------------
def test_bundleArmStreams(ad_blue, ad_red):
    """Test that bundleArmStreams re-bundles the arm streams into one bundle."""
    p = MaroonXSpectrum([ad_blue])
    p.streams['RED'] = [ad_red]
    bundled_list = p.bundleArmStreams(suffix='')

    assert len(bundled_list) == 1, 'Should produce one bundle'

    bundle_ad = bundled_list[0]

    # Verify the bundle structure
    assert len(bundle_ad) == 2, 'Bundle should have 2 extensions'
    assert bundle_ad.filename == ARCHNAME, 'Filename should be restored'

    # Verify both arms are present
    assert 'BLUE' in bundle_ad[0].tags, 'Bundle should contain BLUE arm'
    assert 'RED' in bundle_ad[1].tags, 'Bundle should contain RED arm'
    assert 'BUNDLE' in bundle_ad.tags, 'Bundle should have BUNDLE tag'

    # Verify ORIGNAME is set correctly
    assert bundle_ad.phu.get('ORIGNAME') == ARCHNAME, 'ORIGNAME should be archive name'

    # Verify ARCHNAME was removed (it's now the filename)
    assert 'ARCHNAME' not in bundle_ad.phu, 'ARCHNAME should be removed after bundling'


def test_bundleArmStreams_no_red_stream(ad_blue):
    """Test that bundleArmStreams raises error when RED stream is missing."""
    p = MaroonXSpectrum([ad_blue])

    with pytest.raises(ValueError, match='RED stream not found'):
        p.bundleArmStreams()
