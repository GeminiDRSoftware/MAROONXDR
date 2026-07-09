.. MaroonX DRAGONS documentation master file

==================================
DRAGONS MAROON-X DRP Documentation
==================================

Welcome to the Maroon-X DRP docs. This documentation is intended for
users of the Maroon-X DRP who want to reduce their observations as provided
by the `Gemini Observatory Archive <https://archive.gemini.edu/searchform>`_.

Maroon-X is Gemini North's high-resolution optical echelle spectrograph, 
built for extreme-precision radial-velocity measurements in the search for 
Earth-like planets around M dwarfs. What follows covers the 
DRAGONS-based data reduction pipeline that turns raw Maroon-X frames 
into wavelength-calibrated 1D spectra with barycentric-correction, ready for science.

This documentation is organized in three complementary manuals. 
The **Tutorial** is the primary reference and the fastest way to get going:
it covers installation, calibration database setup, and end-to-end reduction 
workflows. The **User Manual** is a companion reference with in-depth descriptions of 
individual reduction steps (recipes and primitives) and their configuration options.
The **Programmer Manual** documents the pipeline's internals: the DRAGONS instrument
definition, class hierarchy, database integration and tag system. 
It is intended for developers extending or maintaining the DRP.

.. note::

   The **Programmer Manual** is not published online but it comes with the
   repository and can be built locally via ``nox`` sessions when needed. See
   :ref:`maroonx_local_manuals` in the tutorial for details.

Documentation Sections
======================

.. toctree::
   :maxdepth: 2
   :titlesonly:

   tutorials/MAROONXDR_Tutorial/index
   usermanuals/MAROONXDR_UserManual/index




