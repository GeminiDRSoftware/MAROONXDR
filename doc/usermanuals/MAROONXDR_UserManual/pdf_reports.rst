.. pdf_reports.rst

.. _maroonxdr_user_pdf_reports:

**********************
Diagnostic PDF Reports
**********************

Overview
========

Several reduction steps write a diagnostic PDF report to the working
directory alongside their output files. The reports are quality-control
products: each one collects the plots needed to judge whether the step
worked as expected. They are enabled by default through the ``report``
parameter of the primitives listed below.

.. list-table::
   :header-rows: 1
   :widths: 20 25 30 25

   * - Report
     - Primitive
     - Recipe
     - Filename
   * - Background fit
     - :ref:`removeStrayLight <primitive_removeStrayLight>`
     - ``makeProcessedFlat``, ``makeProcessedFlatDFFFF``, ``reduce``
     - ``<input>_backgroundfit.pdf``
   * - Etalon wavelength fit
     - :ref:`fitAndApplyEtalonWls <primitive_fitAndApplyEtalonWls>`
     - ``makeDynamicWavecal``
     - ``<input>_spline.pdf``
   * - Wavelength solution
     - :ref:`applyWavelengthSolution <primitive_applyWavelengthSolution>`
     - ``reduce``
     - ``<input>_spline.pdf``
   * - Fiber combination
     - :ref:`combineFibers <primitive_combineFibers>`
     - ``reduce``
     - ``<input>_fiber6.pdf``
   * - Exposure meter
     - :ref:`barycentricCorrection <primitive_barycentricCorrection>`
     - ``reduce``, ``applyBarycentricCorrection``
     - ``<input>_exposuremeter.pdf``

``<input>`` stands for the name of the FITS file being processed at that
step. Filenames are shown for the default parameters; a few primitives
write variant names for non-default settings, noted in their sections
below.

.. note::

   In a science reduction the background fit report is not produced by
   calling ``removeStrayLight`` directly: the ``reduce`` recipe requests
   straylight removal on the simultaneous calibration fiber through
   ``extractStripes``, which runs ``removeStrayLight`` internally.

Disabling Reports
=================

Each primitive in the table above has its own ``report`` parameter, so
reports are switched off per primitive. On the command line:

.. code-block:: bash

   reduce @flats.list -p removeStrayLight:report=False

In the Python API:

.. code-block:: python

   p.removeStrayLight(report=False)

.. note::

   In a science reduction the background fit report is controlled through
   ``extractStripes``, which forwards the parameter to the straylight
   removal primitive internally:

   .. code-block:: bash

      reduce @science.list -p extractStripes:report=False

Background Fit Report
=====================

Written by :ref:`removeStrayLight <primitive_removeStrayLight>` as
``<input>_backgroundfit.pdf``, for the flat frames in the flat recipes
and for each science frame in ``reduce``.

The seven pages follow the straylight fit through its stages: the
bias-corrected raw frame, the frame with the echelle orders masked, the
first-iteration background model and the corresponding
background-subtracted data, the number of masked pixels per background
mesh cell, the correction-step background model, and the final
background-subtracted frame.

**What to look for.** The order mask should cover every order; an order
missed by the mask leaks into the background model. The background
models should be smooth, and the final background-subtracted frame
should show flat inter-order regions without residual gradients or
outliers.

.. tip::

   The reduction log complements this report. It reports the number of
   pixels rejected during the optimal extraction of each order of the
   science fibers. More than roughly 100 rejected pixels in an order is worth
   inspecting, although high signal-to-noise targets have also been
   observed to produce large rejection counts.

.. todo::
    add example figure

Etalon Wavelength Fit Report
============================

Written by :ref:`fitAndApplyEtalonWls <primitive_fitAndApplyEtalonWls>`
in the ``makeDynamicWavecal`` recipe, as ``<input>_spline.pdf`` named
after the etalon frame (``_spline_symmetrical.pdf`` when
``symmetric_linefits=True``). It shares the ``_spline`` suffix with the
wavelength solution report of the science reduction (next section); the
two are told apart by the frame name (``DEEEE`` versus ``SOOOE``).

The report opens with one etalon dispersion page per etalon-lit fiber,
followed by one residuals page per fiber titled "Etalon residuals after
Etalon-based spline fit", with three panels: residuals against
wavelength, residual RMS per order, and residuals against normalized
detector position.

**What to look for.** The residuals should scatter flat around zero.
Red flags are extremely high residuals and orders whose residual RMS
does not follow the smooth curved trend of the neighboring orders.

.. tip::

   Etalon peak fitting happens in the preceding
   ``getPeaksAndPolynomials`` step, and its failures appear in the
   reduction log as warnings of the form ``Error extracting fiber X
   order Y``. The most common message, an inconsistent number of peak
   minima and maxima, is usually caused by leftover image artifacts or
   cosmic rays.

.. todo::
    add example figure

Wavelength Solution Report
==========================

Written by
:ref:`applyWavelengthSolution <primitive_applyWavelengthSolution>` in
the ``reduce`` recipe, as ``<input>_spline.pdf`` named after the science
frame.

The first page shows the offset of the simultaneous calibration fiber
relative to the etalon frame of the wavelength solution, in three
panels: against wavelength, as the mean fitted offset per order, and
against detector position; this is the drift the correction removes. A second page shows the etalon residuals after the
simultaneous calibration correction, one panel per fiber. The remaining
pages show the residuals after the spline fit for each science fiber and
for the calibration fiber, in the same three-panel layout as the etalon
wavelength fit report.

**What to look for.** The same criteria as the etalon wavelength fit
report: red flags are extremely high residuals and orders that do not
follow the curved trend of the rest.

.. tip::

   The measured drifts are written to the reduction log ("Instrument
   Drift" and "Relative drift measured in Fiber 5", in m/s) and stored
   in the ``INSTRUME_DRIFT`` and ``RELATIVE_DRIFT`` header keywords of
   the reduced file.

.. todo::
    add example figure

Fiber Combination Report
========================

Written by :ref:`combineFibers <primitive_combineFibers>` in the
``reduce`` recipe, as ``<input>_fiber6.pdf`` named after the science
frame (``_fiber7_symmetrical.pdf`` when ``symmetric_linefits=True``).
The name reflects where the result is stored: the combined spectrum is
saved as a virtual fiber 6 (or 7).

The report has one page per echelle order, with four panels: the
intensity of fibers 2, 3 and 4 together with the combined spectrum, the
sigma deviations from the median, the weights (inverse variance), and
the SNR of the combined spectrum.

**What to look for.** The three fibers should reasonably agree; clear
outliers in a single fiber are the red flag.

.. tip::

   For each order the reduction log lists the clipping threshold
   (``kappa_sigma``) and the number of clipped pixels per fiber, and
   warns when the clip count exceeds ``max_clips``. More than 20
   clipped pixels in an order is worth inspecting; the tutorial
   reduction sets ``combineFibers:max_clips=20`` so the log warns at
   that threshold.

.. todo::
    add example figure

Exposure Meter Report
=====================

Written by
:ref:`barycentricCorrection <primitive_barycentricCorrection>` in the
``reduce`` and ``applyBarycentricCorrection`` recipes, as
``<input>_exposuremeter.pdf``.

The single page shows the time series of the two exposure meter
channels (PC and FRD) on dual axes, with the exposure boundaries and
the flux zeropoints of both channels.

**What to look for.** The zeropoints should line up with the actual
flux level outside the exposure, and the flux should be free of
anomalies: no sudden drops or spikes, and especially no drop below the
zeropoint during the exposure.

.. tip::

   The computed barycentric corrections are stored in header keywords
   of the output file: ``BERV_MIDPOINT``, ``BERV_FLUXWEIGHTED_PC`` and
   ``BERV_FLUXWEIGHTED_FRD``, among others. The reduction log records
   the automatic zeropoint determination for both channels and warns
   when it fails, when outliers are replaced, or when gaps are found in
   the exposure meter data.

.. todo::
    add example figure
