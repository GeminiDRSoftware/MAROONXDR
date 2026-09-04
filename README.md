# MAROONXDR

DRAGONS implementation of the data reduction pipeline for the MAROON-X echelle spectrograph at Gemini-North.

[![Testing](https://github.com/GeminiDRSoftware/MAROONXDR/actions/workflows/testing.yml/badge.svg)](https://github.com/GeminiDRSoftware/MAROONXDR/actions/workflows/testing.yml) [![Documentation Status](https://readthedocs.org/projects/maroonxdr/badge/?version=latest)](https://maroonxdr.readthedocs.io/latest/)

> **Status: pre-release / work in progress.** This pipeline is under active
> development. Master darks and synthetic darks, master flats, dynamic
> wavelength calibration, and 2D to 1D echelle
> extraction with drift-corrected wavelength calibration, fiber combination
> and barycentric correction are functional. User and local calibration
> databases are fully implemented.

Full documentation: **https://maroonxdr.readthedocs.io/latest/**

---

## Introduction

MAROON-X is a fiber-fed echelle-spectrograph installed on the Gemini-North Observatory meant to detect Earth-sized planets in the habitable zones of mid- to late-M dwarfs by measuring their radial velocities with 1 m/s radial velocity precision.

The red-optical, fiber-fed echelle-spectrograph has no movable parts and only one mode. At each observed wavelength the 5 fibers are arranged along the cross-dispersion direction. The middle three fibers are pupil-sliced fractions of the on-sky target fiber, the outer two fiber traces are an off-target on-sky fiber and a fiber known as the 'sim. cal. fiber' which can be fed calibration light during science frames (as well as during other frames). All observed light is split with a dichroic to a 'blue arm' (491-670nm) and a 'red arm' (649-920nm). Both arms terminate in CCD detectors with 16-bit 4400x4400 arrays, but the 'blue' detector has 4 reads and the 'red' detector has 2.

## What the pipeline produces

Recipes available in `maroonxdr/maroonx/recipes/sq/`:

| Recipe | Output |
|---|---|
| `recipes_BUNDLE.processBundle` | Split GOA bundle into red/blue arm frames |
| `recipes_DARK.makeProcessedDark` | Master dark |
| `recipes_DARK.makeDarkCoefficients` | Per-pixel dark coefficients fitted from master darks |
| `recipes_DARK.makeSyntheticDarksFromCoeffs` | Synthetic darks at chosen exposure times, built from the dark coefficients |
| `recipes_ECHELLE_SPECT.makeSyntheticDark` | Synthetic dark matching a science exposure |
| `recipes_FLAT_SPECT.makeProcessedFlat` | Master flat with stripe traces and 1D extractions |
| `recipes_FLAT_SPECT.makeBlaze` | Blaze function for each fiber of a master flat |
| `recipes_DYNAMIC_WAVECAL.makeDynamicWavecal` | Etalon wavelength solution used for drift correction (processed wavecal) |
| `recipes_ECHELLE_SPECT.reduce` | Wavelength-calibrated, fiber-combined 1D spectra with barycentric correction |
| `recipes_ECHELLE_SPECT.applyBarycentricCorrection` | Recomputed barycentric correction with target-specific parameters |
| `recipes_ECHELLE_SPECT.exportReducedBundle` | Reduced red and blue arm spectra re-bundled into one file |

Interactive QA variants in `maroonxdr/maroonx/recipes/qa/` (`reduceQA`,
`makeProcessedFlatQA`) display the extracted spectra at a browser-based
Bokeh checkpoint.

## Development installation

The pipeline is currently distributed for development only; there is no
PyPI or conda package yet.

```
git clone https://github.com/GeminiDRSoftware/MAROONXDR.git
cd MAROONXDR
pip install nox
nox -s devenv
source venv/bin/activate
```

`nox -s devenv` creates a Python 3.12 virtualenv, clones DRAGONS into
`./DRAGONS/`, and installs all pipeline + framework dependencies in
editable mode. A `conda` variant is available via `nox -s devconda`.


## Reference lookup files

The pipeline needs a set of static instrument reference files that are not
tracked in git: bad pixel masks (BPM), stripe ID traces (SID), and static
wavelength solutions (WLS). They are distributed with each release as
`lookups_files.zip` on the GitHub releases page.

After installing, extract the archive into the package lookup directory:

```
cd MAROONXDR
unzip lookups_files.zip -d maroonxdr/maroonx/lookups/
```

This places the FITS files under `lookups/BPM/`, `lookups/SID/` and
`lookups/WLS/`, where the pipeline expects them. The development
installation is editable, so no reinstallation is needed.

| Set | Files | Contents |
|---|---|---|
| BPM | `BPM_b_0000.fits`, `BPM_r_0000.fits` | Nominal bad pixel mask per arm |
| SID | `SID_b.fits`, `SID_r.fits` | Nominal stripe ID (fiber and order traces) per arm |
| WLS | `WLSTAT_b.fits`, `WLSTAT_r.fits`, `REFWAVELENGTH_b.fits`, `REFWAVELENGTH_r.fits` | Static wavelength solutions and reference wavelengths per arm |


## Documentation

- **User Manual** — reducing MAROON-X data with this pipeline
- **Tutorial** — step-by-step walkthrough of a full reduction
- **Programmer Manual** — primitives, recipes, internals

The User Manual and the Tutorial are at **https://maroonxdr.readthedocs.io/latest/**.
The Programmer Manual is distributed as a PDF with each release on the
GitHub releases page.

## License

MAROONXDR is distributed under the BSD 3-Clause license, the same license
used by DRAGONS.
