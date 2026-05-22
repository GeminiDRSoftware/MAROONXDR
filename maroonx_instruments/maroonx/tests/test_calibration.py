"""Tests for CalibrationMAROONX calibration association rules.

Tests that a LocalDB can store processed darks and flats, and that
a science-like frame retrieves the correct calibration by arm matching.
"""

import io
import os
import shutil

import astrodata
import pytest

import maroonx_instruments  # noqa - registers AstroDataMAROONX
from recipe_system import cal_service
from recipe_system.config import globalConf

# Re-use the same split files from test_maroonx inputs.
# The create_inputs() below copies them and stamps PROCDARK / PROCFLAT
# to make them appear as processed calibrations.
blue_dark = '20250721T170049Z_DDDDE_b_0300.fits'
red_dark = '20250721T170049Z_DDDDE_r_0300.fits'
blue_flat = '20250701T170101Z_DFFFD_b_0008.fits'
red_flat = '20250701T170101Z_DFFFD_r_0002.fits'
blue_wavecal = '20250717T163124Z_DEEEE_b_0010.fits'

# Processed-stamped copies (created by create_inputs)
proc_blue_dark = 'proc_blue_dark.fits'
proc_red_dark = 'proc_red_dark.fits'
proc_blue_flat = 'proc_blue_flat.fits'
proc_blue_wavecal = 'proc_blue_wavecal.fits'


def _init_caldb(tmp_path):
    """Create and initialize a LocalDB in tmp_path."""
    dbfile = os.path.join(str(tmp_path), 'test.db')
    f = io.StringIO(f"[calibs]\ndatabases = {dbfile} get store")
    globalConf.read_file(f)
    caldb = cal_service.set_local_database()
    caldb.init(wipe=True)
    return caldb


@pytest.mark.preprocessed_data
def test_store_and_retrieve_dark(path_to_inputs, tmp_path):
    """Store a processed blue dark; retrieve it using a blue wavecal as science."""
    caldb = _init_caldb(tmp_path)

    cal_path = os.path.join(path_to_inputs, proc_blue_dark)
    caldb.add_cal(cal_path)

    # Use the raw blue dark as query frame (same arm + same exposure time)
    ad_sci = astrodata.open(os.path.join(path_to_inputs, blue_dark))
    result = caldb.get_calibrations([ad_sci], caltype='processed_dark')
    assert result.files[0] is not None
    assert os.path.basename(result.files[0]) == proc_blue_dark


@pytest.mark.preprocessed_data
def test_store_and_retrieve_flat(path_to_inputs, tmp_path):
    """Store a processed blue flat; retrieve it using a blue wavecal as science."""
    caldb = _init_caldb(tmp_path)

    cal_path = os.path.join(path_to_inputs, proc_blue_flat)
    caldb.add_cal(cal_path)

    ad_sci = astrodata.open(os.path.join(path_to_inputs, blue_wavecal))
    result = caldb.get_calibrations([ad_sci], caltype='processed_flat')
    assert result.files[0] is not None
    assert os.path.basename(result.files[0]) == proc_blue_flat


@pytest.mark.preprocessed_data
def test_no_cross_arm_match(path_to_inputs, tmp_path):
    """A red dark should NOT match a blue science frame."""
    caldb = _init_caldb(tmp_path)

    cal_path = os.path.join(path_to_inputs, proc_red_dark)
    caldb.add_cal(cal_path)

    # Use blue dark as query — same exposure time but wrong arm
    ad_sci = astrodata.open(os.path.join(path_to_inputs, blue_dark))
    result = caldb.get_calibrations([ad_sci], caltype='processed_dark')
    assert result.files[0] is None


@pytest.mark.preprocessed_data
def test_store_and_retrieve_wavecal(path_to_inputs, tmp_path):
    """Store a processed blue wavecal; retrieve it using a blue wavecal as science."""
    caldb = _init_caldb(tmp_path)

    cal_path = os.path.join(path_to_inputs, proc_blue_wavecal)
    caldb.add_cal(cal_path)

    ad_sci = astrodata.open(os.path.join(path_to_inputs, blue_wavecal))
    result = caldb.get_calibrations([ad_sci], caltype='processed_wavecal')
    assert result.files[0] is not None
    assert os.path.basename(result.files[0]) == proc_blue_wavecal


# -- Create inputs -------------------------------------------------------------

def create_inputs():
    """Create processed-stamped calibration files from raw split files.

    Copies raw files from the test_maroonx inputs directory and adds
    PROCDARK or PROCFLAT keywords to make them appear as processed
    calibrations to the FitsStorage ingester.
    """
    source_dir = os.path.join(
        os.environ['DRAGONS_TEST'],
        'maroonx_instruments', 'maroonx', 'test_maroonx', 'inputs',
    )
    dest_dir = os.path.join(
        os.environ['DRAGONS_TEST'],
        'maroonx_instruments', 'maroonx', 'test_calibration', 'inputs',
    )
    os.makedirs(dest_dir, exist_ok=True)

    # Copy raw files needed as query targets
    for name in [blue_wavecal, blue_dark]:
        shutil.copy2(os.path.join(source_dir, name),
                     os.path.join(dest_dir, name))
        print(f'  Copied {name}')

    # Create processed calib files by hacking the header of raw files
    processed_cals = {
        proc_blue_dark:    (blue_dark, 'PROCDARK'),
        proc_red_dark:     (red_dark, 'PROCDARK'),
        proc_blue_flat:    (blue_flat, 'PROCFLAT'),
        proc_blue_wavecal: (blue_wavecal, 'PROCARC'),
    }

    for proc_name, (raw_name, keyword) in processed_cals.items():
        ad = astrodata.open(os.path.join(source_dir, raw_name))
        ad.phu[keyword] = ''
        out_path = os.path.join(dest_dir, proc_name)
        ad.write(out_path, overwrite=True)

        # double check
        ad2 = astrodata.open(out_path)
        assert 'PROCESSED' in ad2.tags, f'{proc_name} missing PROCESSED tag'
        print(f'  Created {proc_name} (tags: {ad2.tags})')


if __name__ == '__main__':
    import sys
    if '--create-inputs' in sys.argv[1:]:
        create_inputs()
    else:
        pytest.main([__file__, '-v'])
