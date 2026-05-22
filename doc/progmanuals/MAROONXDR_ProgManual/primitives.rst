.. _primitives:

**********
Primitives
**********

This section documents all the primitives available in the MAROONX Data Reduction pipeline.
These primitives are organized by functionality and processing stage.

.. note::

   This document provides detailed information about each primitive method,
   including parameters, return values, and usage examples. For a quick
   reference of all primitives, see :ref:`primitives_index`.

2D Image Processing Primitives
===============================

.. note::
   automodule commented out to avoid ReadTheDocs build timeout.

.. .. automodule:: maroonxdr.maroonx.primitives_maroonx_2D
..    :members: stackFlats
..    :no-inherited-members:

Echelle Extraction Primitives
==============================

.. .. automodule:: maroonxdr.maroonx.primitives_maroonx_echelle
..    :members:
..    :no-inherited-members:

Spectrum Processing Primitives
===============================

.. .. automodule:: maroonxdr.maroonx.primitives_maroonx_spectrum
..    :members:
..    :no-inherited-members:

PDF Reports
============

Several primitives produce diagnostic PDF reports to help assess the quality of
each reduction step. Reports are enabled by default via the ``report`` boolean
parameter (default ``True``) on each primitive and are written to the current
working directory. All plotting functions are defined in
``maroonxdr/maroonx/maroonx_plots.py`` and use
``matplotlib.backends.backend_pdf.PdfPages`` for PDF generation.

Primitives That Generate Reports
---------------------------------

``removeStrayLight``
^^^^^^^^^^^^^^^^^^^^

- **Source**: ``primitives_maroonx_2D.py``
- **Output**: ``{input}_backgroundfit.pdf`` (7 pages)
- **Pages**: raw frame, orders masked, background model (1st iteration),
  background subtracted, masked pixels per mesh cell, correction background
  model, final subtracted data
- **Called in**: ``makeProcessedFlat`` recipe
- **Plotting function**: :func:`~maroonxdr.maroonx.maroonx_plots.plot_backgroundfit`

``fitAndApplyEtalonWls``
^^^^^^^^^^^^^^^^^^^^^^^^

- **Source**: ``primitives_maroonx_spectrum.py``
- **Output**: ``spline_[symmetrical_]{input}.pdf``
- **Content**: etalon dispersion plots per fiber, wavelength solution residuals
  after spline fit
- **Called in**: ``makeDynamicWavecal`` recipe
- **Plotting functions**: :meth:`MXEchelleSpectrum.plot_etalon_dispersion`
  (method on spectrum object) and
  :func:`~maroonxdr.maroonx.maroonx_plots.plot_residuals`

``applyWavelengthSolution``
^^^^^^^^^^^^^^^^^^^^^^^^^^^

- **Source**: ``primitives_maroonx_spectrum.py``
- **Output**: ``{input}_spline.pdf``
- **Content**: calibration fiber offset measurements vs wavelength, per-order
  fit residuals, etalon residuals per fiber
- **Called in**: ``reduce`` recipe (``ECHELLE_SPECT``)
- **Plotting functions**:
  :func:`~maroonxdr.maroonx.maroonx_plots.plot_calibfiber_offset`,
  :func:`~maroonxdr.maroonx.maroonx_plots.plot_etalon_residuals`

``combineFibers``
^^^^^^^^^^^^^^^^^

- **Source**: ``primitives_maroonx_spectrum.py``
- **Output**: ``{input}_fiber[6|7_symmetrical].pdf``
- **Content**: one page per echelle order showing intensity comparison
  (fibers 2, 3, 4 + combined), sigma deviations from median, weights
  (1/variance), and SNR
- **Called in**: ``reduce`` recipe (``ECHELLE_SPECT``)
- **Plotting function**:
  :func:`~maroonxdr.maroonx.maroonx_plots.plot_fiber_combination`

``barycentricCorrection``
^^^^^^^^^^^^^^^^^^^^^^^^^

- **Source**: ``primitives_maroonx_spectrum.py``
- **Output**: ``{input}_exposuremeter.pdf`` (1 page)
- **Content**: dual-axis time series of PC and FRD exposure meter channels with
  exposure boundaries and zeropoints
- **Called in**: ``applyBarycentricCorrection`` sub-recipe
- **Plotting function**:
  :func:`~maroonxdr.maroonx.maroonx_plots.plot_exposuremeter`

Plotting Functions Reference
-----------------------------

The following functions in ``maroonxdr/maroonx/maroonx_plots.py`` produce the
figures written into the PDF reports.

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Function
     - Description
   * - ``plot_backgroundfit()``
     - Returns a list of 7 Figures covering each stage of background fitting
       and subtraction.
   * - ``plot_residuals()``
     - Returns 1 Figure with 3 subplots: residuals vs wavelength, residual RMS
       per order, and residuals vs normalized X position.
   * - ``plot_fiber_combination()``
     - Returns 1 Figure with 4 subplots: intensity, sigma deviations from
       median, weights (1/variance), and SNR.
   * - ``plot_calibfiber_offset()``
     - Returns 1 Figure with 3 subplots. Supports overlay mode via the ``fig``
       parameter to plot multiple fibers on the same axes.
   * - ``plot_etalon_residuals()``
     - Returns 1 Figure with 4 subplots. Supports incremental fill via the
       ``plotnumber`` and ``fig`` parameters.
   * - ``plot_exposuremeter()``
     - Returns 1 Figure with dual y-axes (PC and FRD exposure meter channels).

Disabling Reports
------------------

Reports can be disabled by passing ``report=False`` to any primitive. In a
custom recipe this is done directly:

.. code-block:: python

   p.removeStrayLight(report=False)

From the command line, use the DRAGONS parameter override system:

.. code-block:: bash

   reduce @flats.list -p removeStrayLight:report=False

