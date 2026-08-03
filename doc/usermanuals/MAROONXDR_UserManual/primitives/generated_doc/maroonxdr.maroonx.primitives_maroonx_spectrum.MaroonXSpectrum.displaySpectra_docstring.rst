
Display extracted spectra in browser using interactive Bokeh visualization.

This primitive launches an interactive Bokeh server that displays MaroonX
extracted spectra in a web browser. Users can zoom, pan, and inspect
individual orders for quality assessment. The primitive blocks execution
until the user clicks "Submit" or closes the browser window.

The visualizer automatically opens at http://localhost:5006 (or next
available port) and provides interactive plots showing all extracted
orders for selected fibers. If wavelength calibration is available and
show_wavelength=True, spectra are displayed in wavelength space;
otherwise pixel space is used.

Parameters
----------
adinputs : list of :class:`~astrodata.AstroData`
    Input AstroData objects with extracted spectra. Must have
    OPTIMAL_REDUCED_FIBER_* or BOX_REDUCED_FIBER_* extensions.

fibers : list of int
    Fiber numbers to display (e.g., [2, 3, 4] for science fibers,
    [6] for combined fiber, [5] for calibration fiber).
    Default is None.

show_wavelength : bool
    If True and wavelength solution exists (WLS_STATIC_FIBER_* or
    WLS_DYNAMIC_FIBER_* extensions), display spectra vs wavelength (nm).
    If False or no wavelength solution, display vs pixel number.
    Default is False.

Returns
-------
list of :class:`~astrodata.AstroData`
    Unmodified input AstroData objects (this is a visualization-only
    primitive with no data modification).

Notes
-----
The Bokeh server runs in a separate thread and communicates with the
browser via websockets. Browser opens automatically based on user
configuration in ~/.dragons/dragonsrc:

[interactive]
browser = chrome          # or firefox, safari
theme = dark_minimal      # or light_minimal
port_number = 5006        # starting port

If the specified port is occupied, the server automatically tries
subsequent ports up to 65535.

The visualizer provides:
- Dropdown to select fiber
- Dropdown to select individual orders or "All Orders" view
- Dropdown to select extraction type (optimal or box)
- Interactive Bokeh plot with zoom, pan, hover tooltips
- Submit button to continue reduction

Examples
--------
Display science fibers after extraction (pixel space):
>>> p.displaySpectra(fibers=[2, 3, 4])

Display calibrated spectra in wavelength space:
>>> p.displaySpectra(fibers=[2, 3, 4], show_wavelength=True)

Display combined fiber spectrum:
>>> p.displaySpectra(fibers=[6], show_wavelength=True)
