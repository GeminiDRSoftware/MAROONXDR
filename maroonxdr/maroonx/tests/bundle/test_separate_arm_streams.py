"""Unit tests for the separateArmStreams primitive."""

import pytest

from maroonxdr.maroonx.primitives_maroonx_spectrum import MaroonXSpectrum

from . import make_arm


# -- Fixtures ------------------------------------------------------------------
@pytest.fixture()
def ad_arms():
    """A BLUE and a RED arm sharing a single ARCHNAME."""
    archname = 'N00000000M0000.fits'
    return [make_arm('BLUE', archname), make_arm('RED', archname)]


# -- Tests ---------------------------------------------------------------------
def test_separateArmStreams(ad_arms):
    """Test that separateArmStreams sorts the arms into the blue and red streams."""
    p = MaroonXSpectrum(ad_arms)
    blue_list = p.separateArmStreams()

    assert 'RED' in p.streams, 'RED stream should be created'
    red_list = p.streams['RED']

    assert len(blue_list) == 1, 'Should have 1 blue arm'
    assert len(red_list) == 1, 'Should have 1 red arm'

    # Verify separation worked correctly
    assert 'BLUE' in blue_list[0].tags, 'Blue list should have BLUE tag'
    assert 'RED' in red_list[0].tags, 'Red list should have RED tag'


def test_separateArmStreams_no_archname(ad_arms):
    """Test error handling in separateArmStreams when ARCHNAME is missing."""
    # Remove ARCHNAME from one arm, so neither group holds both arms
    ad_arms[0].phu['ARCHNAME'] = None

    p = MaroonXSpectrum(ad_arms)

    with pytest.raises(ValueError, match='No BLUE or RED'):
        p.separateArmStreams()


def test_separateArmStreams_single_arm(ad_arms):
    """Test error handling in separateArmStreams when one arm is passed."""
    # Remove one arm
    adinputs = [ad_arms[0]]

    p = MaroonXSpectrum(adinputs)

    with pytest.raises(ValueError, match='No BLUE or RED'):
        p.separateArmStreams()
