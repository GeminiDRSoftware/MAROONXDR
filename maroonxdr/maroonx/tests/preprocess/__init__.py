"""Staging scripts: download raw data and produce the shared calibrations.

Self-contained factory that fetches raw MaroonX bundles from the Gemini Archive
and reduces them into the shared calibrations (master darks, dark coefficients,
master flats, wavecals, synthetic darks) and reduced science frames that the
regression tests build their inputs from. Run manually, module by module, in
dependency order: bundle -> dark -> flat -> wavecal -> science.

Calibrations are passed explicitly between steps; the scripts do not rely on
caldb association to retrieve them.
"""
