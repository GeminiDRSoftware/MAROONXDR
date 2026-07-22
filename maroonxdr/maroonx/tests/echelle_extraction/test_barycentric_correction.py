"""Unit tests for the barycentricCorrection primitive.

The fast tests exercise the exposure-meter processing helper
(``_exposuremeterStats``) on a synthetic EXPOSUREMETER table: automatic
zeropoint determination, in-exposure selection, outlier replacement, and the
sparse-data fallback. They do not call barycorrpy.

The preprocessed_data test re-runs the full primitive on the blessed reduced
science products and compares the header results against the stored values.
It resolves the target through SIMBAD (as the blessed run did), so it needs
network access; the stored keywords are deleted from the input copy first so
a silent SIMBAD-failure skip cannot pass.

Usage:
    python -m maroonxdr.maroonx.tests.echelle_extraction.test_barycentric_correction --create-inputs
"""

import datetime
import os
from copy import deepcopy

import numpy as np
import pytest
from astropy.table import Table
from astropy.time import Time
from gempy.utils import logutils

import astrodata
import maroonx_instruments  # noqa - registers the MaroonX AstroData class
from maroonxdr.maroonx.primitives_maroonx_spectrum import MaroonXSpectrum

from . import make_echelle_frame

# Synthetic exposure meter: one reading per second, constant baseline outside
# a 60 s exposure and constant flux inside it, so the auto-zeropoint equals
# the baseline and the corrected readings are flat.
UTC_START = Time('2025-01-01T00:10:00', scale='utc')
EXPTIME = 60.0
BASELINE = 100.0
PC_FLUX = 1000.0
FRD_FLUX = 2000.0
# 61 samples fall inside the exposure (both endpoints inclusive)
N_IN_EXPOSURE = 61

# Blessed reduced science products, per arm.
datasets = [
    '20250717T144308Z_SOOOE_b_0300_reduced.fits',
    '20250717T144308Z_SOOOE_r_0300_reduced.fits',
]

# SIMBAD-dependent results: compared as floats, with room for catalog drift
# between the blessed run and the test run. Values are in m/s.
BERV_TOLERANCES = {
    'BERV_MIDPOINT': 1.0,
    'BERV_FLUXWEIGHTED_PC': 1.0,
    'BERV_FLUXWEIGHTED_FRD': 1.0,
    'BERV_DIFFERENCE_PC': 0.1,
    'BERV_DIFFERENCE_FRD': 0.1,
}

# Deterministic results (timing chain and exposure-meter statistics):
# compared exactly as stored.
EXACT_KEYS = [
    'UTC_START', 'UTC_CORRECTION', 'UTC_MIDPOINT',
    'UTC_FLUXWEIGHTED_PC', 'UTC_FLUXWEIGHTED_FRD',
    'JD_UTC_START', 'JD_UTC_MIDPOINT',
    'JD_UTC_FLUXWEIGHTED_PC', 'JD_UTC_FLUXWEIGHTED_FRD',
    'BERV_SIMBAD_TARGET',
    'COUNTS_PC_MIN', 'COUNTS_PC_MAX', 'COUNTS_PC_MEDIAN', 'COUNTS_PC_STD',
    'COUNTS_PC_ZP',
    'COUNTS_FRD_MIN', 'COUNTS_FRD_MAX', 'COUNTS_FRD_MEDIAN', 'COUNTS_FRD_STD',
    'COUNTS_FRD_ZP',
    'SCALEFACTOR',
]


def make_emeter_ad(pc_spike=None, thin_exposure=False):
    """An echelle frame with a synthetic EXPOSUREMETER table.

    Readings run from 6 min before to 6 min after the exposure, one per
    second. With thin_exposure only two samples fall inside the exposure
    window, triggering the sparse-data fallback. pc_spike sets one
    PC reading to the given value.
    """
    ad = make_echelle_frame('RED')

    start = UTC_START.datetime - datetime.timedelta(seconds=360)
    end = UTC_START.datetime + datetime.timedelta(seconds=EXPTIME + 360)
    exposure_end = UTC_START.datetime + datetime.timedelta(seconds=EXPTIME)

    timestamps = []
    pc = []
    frd = []
    t = start
    while t <= end:
        in_exposure = UTC_START.datetime <= t <= exposure_end
        if in_exposure and thin_exposure:
            # Keep only the two endpoint samples inside the exposure
            if t not in (UTC_START.datetime, exposure_end):
                t += datetime.timedelta(seconds=1)
                continue
        timestamps.append(t.isoformat())
        pc.append(PC_FLUX if in_exposure else BASELINE)
        frd.append(FRD_FLUX if in_exposure else BASELINE)
        t += datetime.timedelta(seconds=1)

    pc = np.array(pc)
    if pc_spike is not None:
        # An interior sample (30 s into the exposure)
        spike_iso = (
            UTC_START.datetime + datetime.timedelta(seconds=30)
        ).isoformat()
        pc[timestamps.index(spike_iso)] = pc_spike

    table = Table({
        'Timestamp': timestamps,
        'Flux PC Channel': pc,
        'Flux FRD Channel': np.array(frd),
    })
    table.meta['header'] = {'ZP_PC': BASELINE, 'ZP_FRD': BASELINE}
    ad.EXPOSUREMETER = table
    return ad


# -- Tests ---------------------------------------------------------------------
def test_exposuremeterStats():
    """Auto-zeropoint equals the baseline; corrected readings are flat."""
    ad = make_emeter_ad()
    emeter = MaroonXSpectrum([])._exposuremeterStats(ad, UTC_START, EXPTIME)

    for channel, flux in (('pc', PC_FLUX), ('frd', FRD_FLUX)):
        stats = emeter[channel]['stats']
        assert stats['zeropoint'] == BASELINE
        np.testing.assert_allclose(stats['min'], flux - BASELINE)
        np.testing.assert_allclose(stats['max'], flux - BASELINE)
        np.testing.assert_allclose(stats['median'], flux - BASELINE)
        np.testing.assert_allclose(stats['std'], 0.0, atol=1e-12)

        assert len(emeter[channel]['times']) == N_IN_EXPOSURE
        np.testing.assert_allclose(emeter[channel]['readings'],
                                   flux - BASELINE)


def test_exposuremeterStats_explicit_zeropoints():
    """User-supplied zeropoints bypass the automatic determination."""
    ad = make_emeter_ad()
    emeter = MaroonXSpectrum([])._exposuremeterStats(
        ad, UTC_START, EXPTIME, zp_pc=50.0, zp_frd=25.0
    )

    assert emeter['pc']['stats']['zeropoint'] == 50.0
    assert emeter['frd']['stats']['zeropoint'] == 25.0
    np.testing.assert_allclose(emeter['pc']['readings'], PC_FLUX - 50.0)
    np.testing.assert_allclose(emeter['frd']['readings'], FRD_FLUX - 25.0)


def test_exposuremeterStats_replaces_outlier():
    """A single spiked reading is replaced by the running median."""
    ad = make_emeter_ad(pc_spike=5000.0)
    emeter = MaroonXSpectrum([])._exposuremeterStats(ad, UTC_START, EXPTIME)

    np.testing.assert_allclose(emeter['pc']['readings'], PC_FLUX - BASELINE)


def test_exposuremeterStats_sparse_data_fallback():
    """Two or fewer in-exposure samples abandon the photometric weighting."""
    ad = make_emeter_ad(thin_exposure=True)
    emeter = MaroonXSpectrum([])._exposuremeterStats(ad, UTC_START, EXPTIME)

    for channel in ('pc', 'frd'):
        assert emeter[channel]['times'] is None
        assert np.all(np.isnan(emeter[channel]['readings']))


@pytest.mark.preprocessed_data
@pytest.mark.regression
@pytest.mark.parametrize('filename', datasets)
def test_barycentricCorrection_real_science(
    filename, path_to_inputs, change_working_dir
):
    """BERV results for a real science frame are stable against the ref."""
    ad = astrodata.open(os.path.join(path_to_inputs, filename))

    # Stored output from the reference run
    ref = {key: ad[0].hdr[key] for key in [*BERV_TOLERANCES, *EXACT_KEYS]}

    # Remove the stored results from the input copy: the primitive skips the
    # file on a failed SIMBAD lookup, and the test must fail loudly (missing
    # keywords) in that case instead of comparing the input to itself.
    ad_in = deepcopy(ad)
    for key in ref:
        del ad_in[0].hdr[key]

    with change_working_dir():
        logutils.config(file_name=f'log_{filename.replace(".fits", "")}.txt')
        ad_out = MaroonXSpectrum([]).barycentricCorrection(
            [ad_in], report=False
        )[0]

    for key, atol in BERV_TOLERANCES.items():
        np.testing.assert_allclose(
            float(ad_out[0].hdr[key]), float(ref[key]), atol=atol,
            err_msg=key,
        )
    for key in EXACT_KEYS:
        assert ad_out[0].hdr[key] == ref[key], (
            f'{key}: {ad_out[0].hdr[key]!r} != {ref[key]!r}'
        )


# -- Recipes to create the inputs ----------------------------------------------
def _module_data_dir():
    root = os.environ.get('DRAGONS_TEST')
    if root is None:
        raise RuntimeError('DRAGONS_TEST environment variable not set')
    module_name = os.path.splitext(os.path.basename(__file__))[0]
    return os.path.join(
        root, 'maroonxdr', 'maroonx', 'echelle_extraction', module_name
    )


def create_inputs_recipe():
    """Copy the blessed reduced science products into inputs/."""
    import shutil
    from pathlib import Path

    inputs_dir = Path(_module_data_dir()) / 'inputs'
    inputs_dir.mkdir(parents=True, exist_ok=True)
    preprocessed_dir = Path(os.environ['DRAGONS_TEST']) / 'preprocessed_files'

    for filename in datasets:
        shutil.copy2(preprocessed_dir / filename, inputs_dir / filename)
        print(f'Copied to inputs/:\n    {filename}')


if __name__ == '__main__':
    import sys

    if '--create-inputs' in sys.argv[1:]:
        create_inputs_recipe()
    else:
        pytest.main([__file__])
