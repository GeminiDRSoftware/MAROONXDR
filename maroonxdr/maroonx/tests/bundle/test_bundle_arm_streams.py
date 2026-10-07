"""Unit tests for the bundleArmStreams primitive."""

import logging

import pytest

from maroonxdr.maroonx.primitives_maroonx_spectrum import MaroonXSpectrum

from . import make_arm

# -- Fixtures ------------------------------------------------------------------
ARCHNAME = 'N00000000M0000.fits'

FLAT_FIBER_SETUP = ['Dark', 'Dark', 'Dark', 'Dark', 'Flat lamp']
ETALON_FIBER_SETUP = ['Dark', 'Etalon', 'Etalon', 'Etalon', 'Etalon']
SCIENCE_FIBER_SETUP = ['Sky', 'Target', 'Target', 'Target', 'Etalon']


def bundled_filename(**arm_kwargs):
    """Bundle a blue and red pair built with make_arm; return the output filename."""
    blue = make_arm('BLUE', ARCHNAME, **arm_kwargs)
    red = make_arm('RED', ARCHNAME, **arm_kwargs)
    p = MaroonXSpectrum([blue])
    p.streams['RED'] = [red]
    return p.bundleArmStreams().pop().filename


@pytest.fixture()
def ad_blue():
    """A BLUE frame."""
    return make_arm('BLUE', ARCHNAME)


@pytest.fixture()
def ad_red():
    """A RED frame."""
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


def test_bundleArmStreams_suffix_dark():
    """A processed dark pair derives the _dark suffix."""
    assert bundled_filename() == 'N00000000M0000_dark.fits'


def test_bundleArmStreams_suffix_dark_coefficients():
    """A dark coefficients pair (COEFF_Z0 attached) derives _darkCoefficients."""
    filename = bundled_filename(attach_coeffs=True)
    assert filename == 'N00000000M0000_darkCoefficients.fits'


def test_bundleArmStreams_suffix_synth_dark():
    """A synthetic dark pair derives the _synth_dark suffix."""
    filename = bundled_filename(phu_keywords={'SYNTHETIC_DARK_CREATED': '2025-01-01'})
    assert filename == 'N00000000M0000_synth_dark.fits'


def test_bundleArmStreams_suffix_flat():
    """A processed flat pair derives the _flat suffix."""
    assert bundled_filename(fiber_setup=FLAT_FIBER_SETUP) == 'N00000000M0000_flat.fits'


def test_bundleArmStreams_suffix_arc():
    """A dynamic wavecal (etalon) pair derives the _arc suffix."""
    assert bundled_filename(fiber_setup=ETALON_FIBER_SETUP) == 'N00000000M0000_arc.fits'


def test_bundleArmStreams_suffix_reduced():
    """A reduced science pair derives the _reduced suffix."""
    filename = bundled_filename(fiber_setup=SCIENCE_FIBER_SETUP)
    assert filename == 'N00000000M0000_reduced.fits'


def test_bundleArmStreams_one_pair_per_archname(caplog):
    """Same-arm files sharing an ARCHNAME collapse to one pair; the first stays."""
    caplog.set_level(logging.DEBUG)
    blue_dark = make_arm('BLUE', ARCHNAME)
    red_dark = make_arm('RED', ARCHNAME)
    blue_coeff = make_arm('BLUE', ARCHNAME, attach_coeffs=True)
    red_coeff = make_arm('RED', ARCHNAME, attach_coeffs=True)

    p = MaroonXSpectrum([blue_dark, blue_coeff])
    p.streams['RED'] = [red_dark, red_coeff]
    bundled_list = p.bundleArmStreams()

    assert len(bundled_list) == 1, 'Should produce one bundle per ARCHNAME'
    # The first input pair stays, so the coefficients pair was skipped
    assert bundled_list[0].filename == 'N00000000M0000_dark.fits'
    assert any('Multiple BLUE files for ARCHNAME' in r.message for r in caplog.records)
    assert any('Multiple RED files for ARCHNAME' in r.message for r in caplog.records)


def test_bundleArmStreams_explicit_suffix_overrides():
    """An explicitly passed suffix wins over the tag-derived one."""
    p = MaroonXSpectrum([make_arm('BLUE', ARCHNAME)])
    p.streams['RED'] = [make_arm('RED', ARCHNAME)]
    filename = p.bundleArmStreams(suffix='_custom').pop().filename
    assert filename == 'N00000000M0000_custom.fits'
