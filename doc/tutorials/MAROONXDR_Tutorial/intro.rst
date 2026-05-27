.. intro.rst

.. _maroonx_setup:

************************************
MAROON-X DRP: Setup and Installation
************************************

Introduction
============

This guide walks you through setting up the MAROON-X data reduction pipeline development environment.
The MAROON-X pipeline (``maroonxdr``) is built on top of the DRAGONS framework and provides recipes
and primitives for reducing MAROON-X echelle spectroscopy data.

The setup process uses **nox** (a Python task automation tool) to create a complete development environment with
all necessary dependencies, including the DRAGONS framework, scientific libraries, and the MAROON-X instrument packages.

Prerequisites
=============

Before beginning, ensure you have:

* **Python 3.12** installed on your system
* **Git** for cloning repositories
* **Poetry** for dependency management (will be installed automatically by nox if missing)

.. note:: The MAROON-X pipeline requires Python 3.12 specifically. Using other Python versions may cause compatibility issues.


Installation Steps
==================

Step 1: Clone the MAROONXDR Repository
---------------------------------------

First, clone the MAROON-X DRAGONS repository:

.. code-block:: bash

    cd /path/to/your/projects
    git clone https://github.com/GeminiDRSoftware/MAROONXDR.git
    cd MAROONXDR

Step 2: Install Nox
--------------------

Nox is the task automation tool that manages the development environment setup. Install it using pip:

.. code-block:: bash

    # Install nox globally or in your user directory
    pip install --user nox

    # Or install globally (may require sudo)
    sudo apt install nox

    # Verify installation
    nox --version

Step 3: List Available Nox Sessions
------------------------------------

View all available nox sessions to understand what tasks are available:

.. code-block:: bash

    # List all nox sessions
    nox -l

``noxfile.py`` defines several sessions; the most relevant to a new
developer are:

.. code-block:: text

    - devenv                  -> Create a development environment.
    - devconda                -> Create a conda development environment.
    - download_raws           -> Download MAROON-X raw files for tests.
    - create_inputs           -> Download and create input files for unit tests.
    - unit_tests              -> Run unit tests.
    - regression_tests        -> Run DRAGONS-style regression tests.
    - legacy_regression_tests -> Run legacy pipeline regression tests.
    - docs                    -> Build documentation using Sphinx.

The ``devenv`` session is what we'll use to create the development environment.

Step 4: Create the Development Environment
-------------------------------------------

Run the ``devenv`` nox session to automatically set up the complete development environment:

.. code-block:: bash

    # Create development environment (this may take a few minutes)
    nox -s devenv

.. important:: This command performs the following steps:

    1. Creates a Python 3.12 virtual environment at ``venv/``
    2. Clones the DRAGONS framework into ``DRAGONS/`` directory (if not already present)
    3. Installs DRAGONS in development mode along with FitsStorage from GitHub
    4. Installs all Python dependencies from ``pyproject.toml`` (main, dev, docs, test groups)
    5. Installs ``pytest_dragons`` from GitHub
    6. Installs ``maroonxdr`` and ``maroonx_instruments`` in editable mode
    7. Sets ``DRAGONS_TEST`` environment variable in the activate script

Once the setup completes, activate the virtual environment:

.. code-block:: bash

    source venv/bin/activate

Your shell prompt should now show ``(mx_dev)`` indicating the environment is active.

.. note:: The ``devenv`` session automatically writes ``DRAGONS_TEST`` into the
   activate script, pointing to a ``mx_test/`` directory next to the repository.
   If you need to run legacy regression tests, set ``MAROONX_LEGACY_TEST``
   manually:

    .. code-block:: bash

        export MAROONX_LEGACY_TEST=/path/to/legacy/data

    See :ref:`tests` for the expected directory structure under each variable.


Step 5: Verify Installation
----------------------------

Verify that the installation was successful:

.. code-block:: bash

    # Check that maroonxdr is installed
    python -c "import maroonxdr; print('maroonxdr installed.')"

    # Check that maroonx_instruments is installed
    python -c "import maroonx_instruments; print('maroonx_instruments installed.')"


.. _maroonx_caldb_setup:

Calibration Database Setup
==========================

The DRAGONS calibration database (``caldb``) keeps track of every processed
calibration produced during a reduction and serves it back to later steps
on demand. The reduction tutorials assume ``caldb`` is configured and
initialised before you start. This is a one-time setup per working
directory.

Create the configuration file
-----------------------------

DRAGONS reads its configuration from ``~/.dragons/dragonsrc``. Create the
directory and file:

.. code-block:: bash

    mkdir -p ~/.dragons
    touch ~/.dragons/dragonsrc

Edit ``~/.dragons/dragonsrc`` so it contains a single ``[calibs]`` section
pointing at the calibration database file for the tutorial. We recommend
keeping the database inside the same working directory you will reduce in,
so the tutorial state is self-contained:

.. code-block:: ini

    [calibs]
    databases = /absolute/path/to/science_dir/cal_manager.db get store

The trailing ``get store`` flags tell DRAGONS to **retrieve** calibrations
from this database during reductions and to **store** every newly produced
calibration into it automatically. Both flags are required for the
tutorials to work as written - without ``store``, every processed
calibration would have to be registered manually with ``caldb add``.

.. note:: Use an absolute path for the database file.

Initialise the database
-----------------------

Once the config file is in place, create the (empty) database file with:

.. code-block:: bash

    caldb init

If a database already exists at the configured path, pass ``-w`` 
(``--wipe``) to start over.

Verify the configuration
------------------------

To confirm DRAGONS is reading the file you just wrote and to see which
database it will use, run:

.. code-block:: bash

    caldb config

The output should show ``~/.dragons/dragonsrc`` as the configuration file
and the absolute path of your ``cal_manager.db`` as the configured
database.

For full details on the calibration service, including remote databases
and multi-database setups, see the DRAGONS Recipe System User Manual.


Understanding the Development Environment
==========================================

Directory Structure After Setup
--------------------------------

After running ``nox -s devenv``, your directory structure will look like this:

.. code-block:: text

    MAROONXDR/
    ├── DRAGONS/                      # DRAGONS framework (cloned by nox)
    │   ├── geminidr/                 # Core DRAGONS primitives
    │   ├── gempy/                    # DRAGONS utilities
    │   └── astrodata/                # AstroData framework
    ├── maroonxdr/                    # MAROON-X pipeline implementation
    │   └── maroonx/                  # Primitives, recipes, tests
    ├── maroonx_instruments/          # MAROON-X instrument definitions
    │   └── maroonx/                  # AstroData class, tags, descriptors
    ├── venv/                         # Virtual environment (created by nox)
    ├── doc/                          # Documentation
    ├── noxfile.py                    # Nox configuration
    ├── pyproject.toml                # Package metadata and dependencies
    └── README.md                     # Repository README



Alternative: Conda Environment
===============================

If you prefer using Conda instead of virtualenv, use the ``devconda`` session:

.. code-block:: bash

    # Create conda development environment
    nox -s devconda

    # Activate the conda environment
    conda activate mx_devconda


Next Steps
==========

Once your development environment is set up, you can proceed to:

* :ref:`maroonx_reduction` - Learn the Python API for data reduction
* :ref:`maroonx_reduction_cli` - Learn the command-line interface for data reduction

Both tutorials will guide you through reducing MAROON-X data from raw frames to science-ready spectra.