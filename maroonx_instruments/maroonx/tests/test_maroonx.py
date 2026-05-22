"""Tests for the AstroDataMAROONX class (adclass.py)."""

import os

import astrodata
import pytest
from astrodata import Section

import maroonx_instruments  # noqa - registers AstroDataMAROONX
from maroonx_instruments.maroonx.adclass import AstroDataMAROONX

# -- Test datasets -------------------------------------------------------------
# DARK bundle: N20250721M6125.fits
blue_dark = '20250721T170049Z_DDDDE_b_0300.fits'
red_dark = '20250721T170049Z_DDDDE_r_0300.fits'

# FLAT bundle: N20250701M6126.fits
blue_flat = '20250701T170101Z_DFFFD_b_0008.fits'
red_flat = '20250701T170101Z_DFFFD_r_0002.fits'

# WAVECAL (Etalon) bundle: N20250717M5948.fits
blue_wavecal = '20250717T163124Z_DEEEE_b_0010.fits'
red_wavecal = '20250717T163124Z_DEEEE_r_0004.fits'

# Unsplit bundle (copy of DARK bundle)
bundle_file = 'N20250721M6125.fits'

# Convenience lists
all_blue = [blue_dark, blue_flat, blue_wavecal]
all_red = [red_dark, red_flat, red_wavecal]
all_split = all_blue + all_red


@pytest.mark.parametrize('filename', [blue_dark, red_flat])
def test_is_right_instance(path_to_inputs, filename):
    """astrodata.open returns an AstroDataMAROONX instance."""
    ad = astrodata.open(os.path.join(path_to_inputs, filename))
    assert isinstance(ad, AstroDataMAROONX)


@pytest.mark.parametrize('filename', [blue_dark])
def test_instrument_descriptor(path_to_inputs, filename):
    """instrument() returns 'MAROONX' (hyphen stripped)."""
    ad = astrodata.open(os.path.join(path_to_inputs, filename))
    assert ad.instrument() == 'MAROONX'


@pytest.mark.parametrize('filename', [blue_dark, red_dark])
def test_ad_length_single_arm(path_to_inputs, filename):
    """A split (single-arm) file has exactly one extension."""
    ad = astrodata.open(os.path.join(path_to_inputs, filename))
    assert len(ad) == 1


def test_ad_length_bundle(path_to_inputs):
    """An unsplit bundle file has two extensions."""
    ad = astrodata.open(os.path.join(path_to_inputs, bundle_file))
    assert len(ad) == 2


@pytest.mark.parametrize('filename', all_blue)
def test_tag_blue(path_to_inputs, filename):
    """Split blue files have the BLUE tag."""
    ad = astrodata.open(os.path.join(path_to_inputs, filename))
    assert 'BLUE' in ad.tags and 'RED' not in ad.tags
    

@pytest.mark.parametrize('filename', all_red)
def test_tag_red(path_to_inputs, filename):
    """Split red files have the RED tag."""
    ad = astrodata.open(os.path.join(path_to_inputs, filename))
    assert 'RED' in ad.tags and 'BLUE' not in ad.tags
    

def test_tag_bundle(path_to_inputs):
    """An unsplit bundle has the BUNDLE tag."""
    ad = astrodata.open(os.path.join(path_to_inputs, bundle_file))
    assert 'BUNDLE' in ad.tags and {'BLUE', 'RED'}.isdisjoint(ad.tags)


@pytest.mark.parametrize('filename', [blue_dark, red_dark])
def test_tag_dark(path_to_inputs, filename):
    """Dark files carry DARK and CAL tags."""
    ad = astrodata.open(os.path.join(path_to_inputs, filename))
    assert {'DARK', 'CAL'} <= ad.tags


@pytest.mark.parametrize('filename', [blue_flat, red_flat])
def test_tag_flat(path_to_inputs, filename):
    """Flat files carry FLAT and CAL tags."""
    ad = astrodata.open(os.path.join(path_to_inputs, filename))
    assert {'FLAT', 'CAL'} <= ad.tags


@pytest.mark.parametrize('filename', [blue_wavecal, red_wavecal])
def test_tag_wavecal(path_to_inputs, filename):
    """Wavecal files carry WAVECAL, SPECT, and CAL tags."""
    ad = astrodata.open(os.path.join(path_to_inputs, filename))
    assert {'WAVECAL', 'SPECT', 'CAL'} <= ad.tags


@pytest.mark.parametrize('filename', all_blue)
def test_arm_descriptor_blue(path_to_inputs, filename):
    """arm() returns ['BLUE'] for blue split files."""
    ad = astrodata.open(os.path.join(path_to_inputs, filename))
    assert ad.arm() == ['BLUE']


@pytest.mark.parametrize('filename', all_red)
def test_arm_descriptor_red(path_to_inputs, filename):
    """arm() returns ['RED'] for red split files."""
    ad = astrodata.open(os.path.join(path_to_inputs, filename))
    assert ad.arm() == ['RED']


@pytest.mark.parametrize('filename', [blue_dark, red_dark])
def test_camera_descriptor(path_to_inputs, filename):
    """camera() returns 'BLUE' or 'RED' for split files."""
    ad = astrodata.open(os.path.join(path_to_inputs, filename))
    assert ad.camera() in ('BLUE', 'RED')


@pytest.mark.parametrize('filename', [blue_dark, red_dark])
def test_fiber_setup_returns_list(path_to_inputs, filename):
    """fiber_setup() returns a 5-element list."""
    ad = astrodata.open(os.path.join(path_to_inputs, filename))
    fs = ad.fiber_setup()
    assert isinstance(fs, list)
    assert len(fs) == 5
    assert fs == ['Dark', 'Dark', 'Dark', 'Dark', 'Etalon']


@pytest.mark.parametrize('filename', [blue_dark, red_dark])
def test_fiber_setup_short(path_to_inputs, filename):
    """fiber_setup(short=True) returns a 5-char string."""
    ad = astrodata.open(os.path.join(path_to_inputs, filename))
    short = ad.fiber_setup(short=True)
    assert isinstance(short, str)
    assert len(short) == 5
    assert short == 'DDDDE'


@pytest.mark.parametrize('filename', [blue_dark])
def test_array_name_blue(path_to_inputs, filename):
    """array_name() returns the 4 blue quadrant names."""
    ad = astrodata.open(os.path.join(path_to_inputs, filename))
    assert ad.array_name() == [['Q1', 'Q2', 'Q3', 'Q4']]


@pytest.mark.parametrize('filename', [red_dark])
def test_array_name_red(path_to_inputs, filename):
    """array_name() returns the 2 red detector names."""
    ad = astrodata.open(os.path.join(path_to_inputs, filename))
    assert ad.array_name() == [['R1', 'R2']]


@pytest.mark.parametrize('filename', [blue_flat, bundle_file])
def test_gain(path_to_inputs, filename):
    """gain() returns a nested list of floats."""
    ad = astrodata.open(os.path.join(path_to_inputs, filename))
    gains = ad.gain()
    assert isinstance(gains, list)
    assert len(gains) == len(ad)
    assert all(isinstance(v, float) for v in gains[0])


@pytest.mark.parametrize('filename', [red_dark, bundle_file])
def test_read_noise(path_to_inputs, filename):
    """read_noise() returns a nested list of floats."""
    ad = astrodata.open(os.path.join(path_to_inputs, filename))
    rn = ad.read_noise()
    assert isinstance(rn, list)
    assert len(rn) == len(ad)
    assert all(isinstance(v, float) for v in rn[0])


@pytest.mark.parametrize('filename', [blue_dark, red_dark])
def test_data_section_returns_sections(path_to_inputs, filename):
    """data_section() returns Section objects."""
    ad = astrodata.open(os.path.join(path_to_inputs, filename))
    ds = ad.data_section()
    assert isinstance(ds, list)
    assert all(isinstance(s, Section) for s in ds[0])


@pytest.mark.parametrize('filename', [blue_dark, red_dark])
def test_exposure_time(path_to_inputs, filename):
    """exposure_time() returns a numeric value."""
    ad = astrodata.open(os.path.join(path_to_inputs, filename))
    assert isinstance(ad.exposure_time(), (int, float))


# -- Create inputs -------------------------------------------------------------
bundles_needed = [
    'N20250721M6125.fits',   # DARK
    'N20250701M6126.fits',   # FLAT
    'N20250717M5948.fits',   # WAVECAL
]


def create_inputs():
    """Split raw bundles into single-arm files for the AstroData tests.

    Reads raw bundles from $DRAGONS_TEST/raw_files/ (populated by the
    download_raws nox session) and writes split-arm files into the
    test_maroonx inputs directory.

    Run with:
        python -m maroonx_instruments.maroonx.tests.test_maroonx --create-inputs
    """
    from maroonxdr.maroonx.primitives_maroonx_2D import MAROONX

    raw_dir = os.path.join(os.environ['DRAGONS_TEST'], 'raw_files')
    input_path = os.path.join(
        os.environ['DRAGONS_TEST'],
        'maroonx_instruments', 'maroonx', 'test_maroonx', 'inputs',
    )
    os.makedirs(input_path, exist_ok=True)

    for filename in bundles_needed:
        raw_path = os.path.join(raw_dir, filename)
        if not os.path.isfile(raw_path):
            print(f'  Skipping {filename}: not in {raw_dir}')
            continue

        ad = astrodata.open(raw_path)

        # Keep one unsplit bundle for bundle tests
        if filename == bundles_needed[0]:
            out = os.path.join(input_path, ad.filename)
            ad.write(out, overwrite=True)
            print(f'  Wrote bundle {ad.filename}')

        p = MAROONX([ad])
        split_ads = p.splitBundle()
        for ad_arm in split_ads:
            ad_arm.write(
                os.path.join(input_path, ad_arm.filename), overwrite=True
            )
            print(f'  Wrote {ad_arm.filename}')


if __name__ == '__main__':
    import sys
    if '--create-inputs' in sys.argv[1:]:
        create_inputs()
    else:
        pytest.main([__file__, '-v'])
