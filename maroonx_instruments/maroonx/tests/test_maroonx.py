"""Tests for the AstroDataMAROONX class (adclass.py).

All input frames are synthetic: minimal headers, tiny data arrays. The
fixture writes them to a temporary directory and the tests reopen them with
astrodata.open.
"""

import os

import astrodata
import pytest
from astrodata import Section

import maroonx_instruments  # noqa - registers AstroDataMAROONX
from maroonx_instruments.maroonx.adclass import AstroDataMAROONX
from . import (
    DARK_FIBER_SETUP,
    ETALON_FIBER_SETUP,
    FLAT_FIBER_SETUP,
    make_bundle,
    make_frame,
)

# -- Test datasets -------------------------------------------------------------
# Filenames produced by make_frame / make_bundle for each (arm, setup) pair.
blue_dark = '00000000T000000Z_DDDDE_b_0300.fits'
red_dark = '00000000T000000Z_DDDDE_r_0300.fits'

blue_flat = '00000000T000000Z_DFFFD_b_0300.fits'
red_flat = '00000000T000000Z_DFFFD_r_0300.fits'

blue_wavecal = '00000000T000000Z_DEEEE_b_0300.fits'
red_wavecal = '00000000T000000Z_DEEEE_r_0300.fits'

bundle_file = 'N00000000M0000.fits'

# Convenience lists
all_blue = [blue_dark, blue_flat, blue_wavecal]
all_red = [red_dark, red_flat, red_wavecal]
all_split = all_blue + all_red


@pytest.fixture(scope='module')
def inputs(tmp_path_factory):
    """Write the synthetic frames to disk, once per module."""
    path = tmp_path_factory.mktemp('test_maroonx_inputs')

    setups = [DARK_FIBER_SETUP, FLAT_FIBER_SETUP, ETALON_FIBER_SETUP]
    for arm in ('BLUE', 'RED'):
        for fiber_setup in setups:
            ad = make_frame(arm, fiber_setup=fiber_setup)
            ad.write(str(path / ad.filename))

    ad = make_bundle()
    ad.write(str(path / ad.filename))

    return str(path)


@pytest.mark.parametrize('filename', [blue_dark, red_flat])
def test_is_right_instance(inputs, filename):
    """astrodata.open returns an AstroDataMAROONX instance."""
    ad = astrodata.open(os.path.join(inputs, filename))
    assert isinstance(ad, AstroDataMAROONX)


@pytest.mark.parametrize('filename', [blue_dark])
def test_instrument_descriptor(inputs, filename):
    """instrument() returns 'MAROONX' (hyphen stripped)."""
    ad = astrodata.open(os.path.join(inputs, filename))
    assert ad.instrument() == 'MAROONX'


@pytest.mark.parametrize('filename', [blue_dark, red_dark])
def test_ad_length_single_arm(inputs, filename):
    """A split (single-arm) file has exactly one extension."""
    ad = astrodata.open(os.path.join(inputs, filename))
    assert len(ad) == 1


def test_ad_length_bundle(inputs):
    """An unsplit bundle file has two extensions."""
    ad = astrodata.open(os.path.join(inputs, bundle_file))
    assert len(ad) == 2


@pytest.mark.parametrize('filename', all_blue)
def test_tag_blue(inputs, filename):
    """Split blue files have the BLUE tag."""
    ad = astrodata.open(os.path.join(inputs, filename))
    assert 'BLUE' in ad.tags and 'RED' not in ad.tags


@pytest.mark.parametrize('filename', all_red)
def test_tag_red(inputs, filename):
    """Split red files have the RED tag."""
    ad = astrodata.open(os.path.join(inputs, filename))
    assert 'RED' in ad.tags and 'BLUE' not in ad.tags


def test_tag_bundle(inputs):
    """An unsplit bundle has the BUNDLE tag."""
    ad = astrodata.open(os.path.join(inputs, bundle_file))
    assert 'BUNDLE' in ad.tags and {'BLUE', 'RED'}.isdisjoint(ad.tags)


@pytest.mark.parametrize('filename', [blue_dark, red_dark])
def test_tag_dark(inputs, filename):
    """Dark files carry DARK and CAL tags."""
    ad = astrodata.open(os.path.join(inputs, filename))
    assert {'DARK', 'CAL'} <= ad.tags


@pytest.mark.parametrize('filename', [blue_flat, red_flat])
def test_tag_flat(inputs, filename):
    """Flat files carry FLAT and CAL tags."""
    ad = astrodata.open(os.path.join(inputs, filename))
    assert {'FLAT', 'CAL'} <= ad.tags


@pytest.mark.parametrize('filename', [blue_wavecal, red_wavecal])
def test_tag_wavecal(inputs, filename):
    """Wavecal files carry WAVECAL, SPECT, and CAL tags."""
    ad = astrodata.open(os.path.join(inputs, filename))
    assert {'WAVECAL', 'SPECT', 'CAL'} <= ad.tags


@pytest.mark.parametrize('filename', all_blue)
def test_arm_descriptor_blue(inputs, filename):
    """arm() returns ['BLUE'] for blue split files."""
    ad = astrodata.open(os.path.join(inputs, filename))
    assert ad.arm() == ['BLUE']


@pytest.mark.parametrize('filename', all_red)
def test_arm_descriptor_red(inputs, filename):
    """arm() returns ['RED'] for red split files."""
    ad = astrodata.open(os.path.join(inputs, filename))
    assert ad.arm() == ['RED']


@pytest.mark.parametrize('filename', [blue_dark, red_dark])
def test_camera_descriptor(inputs, filename):
    """camera() returns 'BLUE' or 'RED' for split files."""
    ad = astrodata.open(os.path.join(inputs, filename))
    assert ad.camera() in ('BLUE', 'RED')


@pytest.mark.parametrize('filename', [blue_dark, red_dark])
def test_fiber_setup_returns_list(inputs, filename):
    """fiber_setup() returns a 5-element list."""
    ad = astrodata.open(os.path.join(inputs, filename))
    fs = ad.fiber_setup()
    assert isinstance(fs, list)
    assert len(fs) == 5
    assert fs == ['Dark', 'Dark', 'Dark', 'Dark', 'Etalon']


@pytest.mark.parametrize('filename', [blue_dark, red_dark])
def test_fiber_setup_short(inputs, filename):
    """fiber_setup(short=True) returns a 5-char string."""
    ad = astrodata.open(os.path.join(inputs, filename))
    short = ad.fiber_setup(short=True)
    assert isinstance(short, str)
    assert len(short) == 5
    assert short == 'DDDDE'


@pytest.mark.parametrize('filename', [blue_dark])
def test_array_name_blue(inputs, filename):
    """array_name() returns the 4 blue quadrant names."""
    ad = astrodata.open(os.path.join(inputs, filename))
    assert ad.array_name() == [['Q1', 'Q2', 'Q3', 'Q4']]


@pytest.mark.parametrize('filename', [red_dark])
def test_array_name_red(inputs, filename):
    """array_name() returns the 2 red detector names."""
    ad = astrodata.open(os.path.join(inputs, filename))
    assert ad.array_name() == [['R1', 'R2']]


@pytest.mark.parametrize('filename', [blue_flat, bundle_file])
def test_gain(inputs, filename):
    """gain() returns a nested list of floats."""
    ad = astrodata.open(os.path.join(inputs, filename))
    gains = ad.gain()
    assert isinstance(gains, list)
    assert len(gains) == len(ad)
    assert all(isinstance(v, float) for v in gains[0])


@pytest.mark.parametrize('filename', [red_dark, bundle_file])
def test_read_noise(inputs, filename):
    """read_noise() returns a nested list of floats."""
    ad = astrodata.open(os.path.join(inputs, filename))
    rn = ad.read_noise()
    assert isinstance(rn, list)
    assert len(rn) == len(ad)
    assert all(isinstance(v, float) for v in rn[0])


@pytest.mark.parametrize('filename', [blue_dark, red_dark])
def test_data_section_returns_sections(inputs, filename):
    """data_section() returns Section objects."""
    ad = astrodata.open(os.path.join(inputs, filename))
    ds = ad.data_section()
    assert isinstance(ds, list)
    assert all(isinstance(s, Section) for s in ds[0])


@pytest.mark.parametrize('filename', [blue_dark, red_dark])
def test_exposure_time(inputs, filename):
    """exposure_time() returns a numeric value."""
    ad = astrodata.open(os.path.join(inputs, filename))
    assert isinstance(ad.exposure_time(), (int, float))
