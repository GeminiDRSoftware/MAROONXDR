"""
Interactive Bokeh visualizer for MaroonX extracted spectra.

This module provides a web-based spectrum viewer that displays MaroonX
echelle spectra with interactive zoom, pan, and order selection capabilities.
"""

import numpy as np
from bokeh.plotting import figure
from bokeh.models import ColumnDataSource, Select, Div, Button, Range1d, Tabs, TabPanel
from bokeh.layouts import column, row

from geminidr.interactive.interactive import PrimitiveVisualizer, UIParameters


class SpectrumVisualizer(PrimitiveVisualizer):
    """
    Interactive Bokeh visualizer for MaroonX extracted spectra.

    This visualizer displays extracted echelle spectra in a web browser,
    allowing users to:
    - Select which file to view (if multiple files)
    - Select which arm to view (Blue/Red for bundles)
    - Select which fiber to view
    - Select individual orders or view all orders stacked
    - Zoom and pan interactively
    - Inspect wavelength/flux values with hover tooltips

    The visualizer automatically detects whether wavelength calibration is
    available and adjusts the x-axis accordingly (wavelength vs pixel).

    Parameters
    ----------
    all_files_data : list of dict
        List of dictionaries, one per file. Each dictionary contains:
        - 'filename': str - Name of the file
        - 'arms_data': dict - Dictionary mapping arm names to spectra data
            For single arm: {'Data': {fiber: {...}}}
            For multiple arms: {'Blue': {fiber: {...}}, 'Red': {fiber: {...}}}
            Each fiber entry contains:
            - 'optimal_flux': 2D numpy array (n_orders × n_pixels) or None
            - 'box_flux': 2D numpy array (n_orders × n_pixels) or None
            - 'orders': 1D array of order numbers
            - 'wavelength': 2D array (n_orders × n_pixels) or None
    title : str
        Title for the visualizer window
    primitive_name : str
        Name of the calling primitive (for display)
    """

    def __init__(
        self,
        all_files_data,
        title="MaroonX Spectra",
        primitive_name="displaySpectra",
        **kwargs
    ):
        # Create empty UIParameters if not provided (read-only visualizer)
        if 'ui_params' not in kwargs:
            kwargs['ui_params'] = UIParameters()

        # Set filename_info based on number of files
        if len(all_files_data) == 1:
            filename_info = all_files_data[0]['filename']
        else:
            filename_info = f"{len(all_files_data)} files"

        super().__init__(
            title=title,
            primitive_name=primitive_name,
            filename_info=filename_info,
            **kwargs
        )

        self.all_files_data = all_files_data
        self.user_satisfied = True  # Always satisfied (read-only viewer)

        # UI elements
        self.fiber_select = None
        self.order_select = None
        self.extraction_select = None
        self.plot = None
        self.data_source = None

    def visualize(self, doc):
        """
        Build and display the Bokeh UI.

        Parameters
        ----------
        doc : bokeh.document.Document
            Bokeh document for the web page
        """
        if len(self.all_files_data) > 1:
            # Multiple files - create file tabs at top level
            file_tabs = []
            for file_data in self.all_files_data:
                filename = file_data['filename']
                arms_data = file_data['arms_data']
                file_content = self._create_file_panel(arms_data, filename)
                file_tab = TabPanel(child=file_content, title=filename)
                file_tabs.append(file_tab)

            tabs_widget = Tabs(tabs=file_tabs)
            layout = column(tabs_widget, self.submit_button)
            doc.add_root(layout)
        else:
            # Single file - create arm tabs directly
            file_data = self.all_files_data[0]
            arms_data = file_data['arms_data']
            filename = file_data['filename']
            layout = self._create_file_panel(arms_data, filename)
            doc.add_root(column(layout, self.submit_button))

    def _create_file_panel(self, arms_data, filename):
        """
        Create a panel for a single file with arm tabs if needed.

        Parameters
        ----------
        arms_data : dict
            Dictionary mapping arm names to spectra data
        filename : str
            Name of the file

        Returns
        -------
        bokeh.layouts.column or Tabs
            Layout containing the arm panels
        """
        # Check if we have multiple arms (create tabs) or single arm
        if len(arms_data) > 1:
            # Multiple arms - create arm tabs
            arm_tabs = []
            for arm_name in sorted(arms_data.keys()):
                spectra_data = arms_data[arm_name]
                tab_content = self._create_arm_panel(spectra_data, arm_name, filename)
                tab = TabPanel(child=tab_content, title=arm_name)
                arm_tabs.append(tab)

            return Tabs(tabs=arm_tabs)
        else:
            # Single arm - no tabs
            arm_name = list(arms_data.keys())[0]
            spectra_data = arms_data[arm_name]
            return self._create_arm_panel(spectra_data, arm_name, filename)

    def _create_arm_panel(self, spectra_data, arm_name, filename=None):
        """
        Create a panel for a single arm with fiber/order/extraction controls.

        Parameters
        ----------
        spectra_data : dict
            Dictionary mapping fiber numbers to spectrum data
        arm_name : str
            Name of the arm (e.g., 'Blue', 'Red', 'Data')
        filename : str, optional
            Name of the file being displayed

        Returns
        -------
        bokeh.layouts.column
            Layout containing the controls and plot
        """
        # Create fiber selection dropdown
        fiber_options = [(str(f), f"Fiber {f}") for f in sorted(spectra_data.keys())]
        initial_fiber = sorted(spectra_data.keys())[0]

        fiber_select = Select(
            title="Fiber:",
            value=str(initial_fiber),
            options=fiber_options,
            width=150
        )

        # Create order selection dropdown
        orders = spectra_data[initial_fiber].get('orders', [])
        order_options = [("all", "All Orders (stacked)")] + [
            (str(int(o)), f"Order {int(o)}") for o in orders
        ]

        order_select = Select(
            title="Order:",
            value="all",
            options=order_options,
            width=200
        )

        # Create extraction type selection dropdown
        # Determine which extraction types are available across all fibers
        available_extractions = set()
        for fib_data in spectra_data.values():
            if 'optimal_flux' in fib_data and fib_data['optimal_flux'] is not None:
                available_extractions.add('optimal')
            if 'box_flux' in fib_data and fib_data['box_flux'] is not None:
                available_extractions.add('box')

        extraction_options = []
        if 'optimal' in available_extractions:
            extraction_options.append(('optimal', 'Optimal'))
        if 'box' in available_extractions:
            extraction_options.append(('box', 'Box'))

        # Default to optimal if available, otherwise box
        initial_extraction = 'optimal' if 'optimal' in available_extractions else 'box'

        extraction_select = Select(
            title="Extraction:",
            value=initial_extraction,
            options=extraction_options,
            width=150
        )

        # Create info panel
        fiber_list = ', '.join([str(f) for f in sorted(spectra_data.keys())])
        extraction_info = []
        if 'optimal' in available_extractions:
            extraction_info.append("<span style='color:navy'>■ Optimal</span>")
        if 'box' in available_extractions:
            extraction_info.append("<span style='color:darkorange'>■ Box</span>")

        # Check if this arm has wavelength data
        has_wavelength = any([spectra_data[f].get('wavelength') is not None for f in spectra_data.keys()])

        # Build info text with optional filename
        info_lines = []
        if filename:
            info_lines.append(f"<b>File:</b> {filename}")
        info_lines.extend([
            f"<b>Arm:</b> {arm_name}",
            f"<b>Wavelength:</b> {'Calibrated' if has_wavelength else 'Pixel space'}",
            f"<b>Fibers:</b> {fiber_list}",
            f"<b>Available:</b> {', '.join(extraction_info)}"
        ])
        info_text = "<br>".join(info_lines)
        info_div = Div(text=info_text, width=400, height=120)

        # Create the plot
        x_label = "Wavelength (nm)" if has_wavelength else "Pixel"
        plot = figure(
            width=1000,
            height=600,
            title=f"{arm_name} - Fiber {initial_fiber}",
            x_axis_label=x_label,
            y_axis_label="Flux (counts)",
            tools="pan,box_zoom,wheel_zoom,reset,save",
            active_scroll="wheel_zoom",
        )

        # Initialize data source
        data_source = ColumnDataSource(data=dict(x=[], y=[]))
        # Store the line renderer so we can update its color
        line_renderer = plot.line(
            'x', 'y', source=data_source, line_width=1.5, color='navy'
        )

        # Update plot with initial data
        self._update_plot_data(
            spectra_data, plot, data_source, line_renderer,
            initial_fiber, "all", initial_extraction
        )

        # Add callbacks with closures to capture local variables
        def on_fiber_change(attr, old, new):
            fiber = int(new)
            orders = spectra_data[fiber].get('orders', [])
            order_options = [("all", "All Orders (stacked)")] + [
                (str(int(o)), f"Order {int(o)}") for o in orders
            ]
            order_select.options = order_options
            order_select.value = "all"
            extraction = extraction_select.value
            self._update_plot_data(
                spectra_data, plot, data_source, line_renderer,
                fiber, "all", extraction
            )

        def on_order_change(attr, old, new):
            fiber = int(fiber_select.value)
            extraction = extraction_select.value
            self._update_plot_data(
                spectra_data, plot, data_source, line_renderer,
                fiber, new, extraction
            )

        def on_extraction_change(attr, old, new):
            fiber = int(fiber_select.value)
            order = order_select.value
            self._update_plot_data(
                spectra_data, plot, data_source, line_renderer,
                fiber, order, new
            )

        fiber_select.on_change('value', on_fiber_change)
        order_select.on_change('value', on_order_change)
        extraction_select.on_change('value', on_extraction_change)

        # Create layout
        controls = row(fiber_select, order_select, extraction_select, info_div)
        return column(controls, plot)

    def _update_plot_data(self, spectra_data, plot, data_source, line_renderer, fiber, order_str, extraction_type):
        """
        Update plot with new fiber/order/extraction data.

        Parameters
        ----------
        spectra_data : dict
            Dictionary mapping fiber numbers to spectrum data
        plot : bokeh.plotting.figure
            The plot object to update
        data_source : ColumnDataSource
            The data source for the plot
        line_renderer : GlyphRenderer
            The line glyph renderer
        fiber : int
            Fiber number to display
        order_str : str
            Order to display ('all' for stacked view, or specific order number)
        extraction_type : str
            Extraction type to display ('optimal' or 'box')
        """
        fiber_data = spectra_data[fiber]
        orders = fiber_data.get('orders')
        wavelength = fiber_data.get('wavelength')

        # Get the appropriate flux array based on extraction type
        if extraction_type == 'optimal':
            flux = fiber_data.get('optimal_flux')
            if flux is None:
                # Fall back to box if optimal not available
                flux = fiber_data.get('box_flux')
                extraction_type = 'box'
        else:
            flux = fiber_data.get('box_flux')
            if flux is None:
                # Fall back to optimal if box not available
                flux = fiber_data.get('optimal_flux')
                extraction_type = 'optimal'

        # Set line color based on extraction type
        # Optimal: navy blue, Box: orange
        line_color = 'navy' if extraction_type == 'optimal' else 'darkorange'
        line_renderer.glyph.line_color = line_color

        if order_str == "all":
            # Stacked view of all orders
            x_data = []
            y_data = []

            # Stack orders with vertical offset
            for i, order in enumerate(orders):
                order_idx = i
                spectrum = flux[order_idx]

                # Get x-axis values
                if wavelength is not None:
                    x = wavelength[order_idx]
                else:
                    x = np.arange(len(spectrum))

                # Normalize and offset for stacking
                median_flux = np.nanmedian(spectrum)
                if median_flux > 0:
                    norm_spectrum = spectrum / median_flux
                else:
                    norm_spectrum = spectrum

                offset = i * 1.5  # Vertical spacing between orders
                y = norm_spectrum + offset

                # Filter out NaNs
                valid = np.isfinite(x) & np.isfinite(y)
                x_data.extend(x[valid])
                y_data.extend(y[valid])

                # Add separator between orders
                if i < len(orders) - 1:
                    x_data.append(np.nan)
                    y_data.append(np.nan)

            extr_label = extraction_type.capitalize()
            plot.title.text = f"All Orders - Fiber {fiber} ({extr_label}, stacked, normalized)"
            plot.yaxis.axis_label = "Normalized + Offset"
            plot.xaxis.axis_label = "Wavelength (nm)" if wavelength is not None else "Pixel"

        else:
            # Single order view
            order = int(order_str)
            order_idx = orders.tolist().index(order)
            spectrum = flux[order_idx]

            # Get x-axis values
            if wavelength is not None:
                x_data = wavelength[order_idx]
            else:
                x_data = np.arange(len(spectrum))

            y_data = spectrum

            # Filter out NaNs
            valid = np.isfinite(x_data) & np.isfinite(y_data)
            x_data = x_data[valid]
            y_data = y_data[valid]

            extr_label = extraction_type.capitalize()
            plot.title.text = f"Order {order} - Fiber {fiber} ({extr_label})"
            plot.yaxis.axis_label = "Flux (counts)"
            plot.xaxis.axis_label = "Wavelength (nm)" if wavelength is not None else "Pixel"

        # Update data source
        data_source.data = dict(x=x_data, y=y_data)

        # Auto-scale axes - replace the range objects entirely to force update
        x_min, x_max = np.nanmin(x_data), np.nanmax(x_data)
        x_range = x_max - x_min
        x_start = x_min - 0.02 * x_range
        x_end = x_max + 0.02 * x_range
        plot.x_range = Range1d(start=x_start, end=x_end)

        y_min, y_max = np.nanmin(y_data), np.nanmax(y_data)
        y_range = y_max - y_min
        y_start = y_min - 0.02 * y_range
        y_end = y_max + 0.02 * y_range
        plot.y_range = Range1d(start=y_start, end=y_end)

    def submit_button_handler(self):
        """
        Handle submit button click.

        For a read-only visualizer, simply accept and close.
        """
        self.user_satisfied = True
        self.submit_button.disabled = True

    def results(self):
        """Return results (None for read-only viewer)."""
        return None
