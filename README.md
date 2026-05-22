# MAROONXDR

DRAGONS implementation of the data reduction pipeline for the MAROON-X echelle spectrograph at Gemini-North.

> **Status: pre-release / work in progress.** This pipeline is under active
> development. Master darks, master flats, dynamic wavelength calibration,
> and 2D to 1D echelle extraction with barycentric correction are functional.
> Full CalDB integration is still in progress.

Full documentation: **https://maroonxdr.readthedocs.io/latest/**

---

## Introduction

MAROON-X is a fiber-fed echelle-spectrograph installed on the Gemini-North Observatory meant to detect Earth-sized planets in the habitable zones of mid- to late-M dwarves by measuring their radial velocities with 1 m/s radial velocity precision.

The red-optical, fiber-fed echelle-spectrograph has no movable parts and only one mode. At each observed wavelength the 5 fibers are arranged along the cross-dispersion direction. The middle three fibers are pupil-sliced fractions of the on-sky target fiber, the outer two fiber traces are an off-target on-sky fiber and a fiber known as the 'sim. cal. fiber' which can be fed calibration light during science frames (as well as during other frames). All observed light is split with a dichroic to a ‘blue arm’ (491-670nm) and a ‘red arm’ (649-920nm). Both arms terminate in CCD detectors with 16-bit 4400x4400 arrays, but the 'blue' detector has 4 reads and the 'red' detector has 2.

## What the pipeline produces

Recipes available in `maroonxdr/maroonx/recipes/sq/`:

| Recipe | Output |
|---|---|
| `recipes_BUNDLE.processBundle` | Split GOA bundle into red/blue arm frames |
| `recipes_DARK.makeProcessedDark` | Master dark |
| `recipes_FLAT_SPECT.makeProcessedFlat` | Master flat with stripe traces |
| `recipes_DYNAMIC_WAVECAL.makeDynamicWavecal` | Per-night etalon wavelength solution |
| `recipes_ECHELLE_SPECT.reduce` | Optimally extracted 1D spectra |
| `recipes_ECHELLE_SPECT.applyBarycentricCorrection` | Compute barycentric-correction |



## Development installation

The pipeline is currently distributed for development only — there is no
PyPI/conda package yet.

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


## Documentation

- **User Manual** — reducing MAROON-X data with this pipeline
- **Tutorial** — step-by-step walkthrough of a full reduction
- **Programmer Manual** — primitives, recipes, internals

All three are at **https://maroonxdr.readthedocs.io/latest/**.
