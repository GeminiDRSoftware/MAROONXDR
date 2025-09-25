"""Noxfile for the MAROONXDR package.

For a list of sessions and their descriptions, run:

.. code-block::
    nox -l

Please note that, to prevent undesirable execution, there are no default
sessions, so running ``nox`` in isolation will do nothing.
"""

import re
from pathlib import Path
import os

import nox

nox.options.sessions = []
# nox.options.default_venv_backend = 'conda'
nox.options.error_on_external_run = True

# Dragons installation resources
DRAGONS_URL = R'https://github.com/GeminiDRSoftware/DRAGONS'
CALMGR_URL = R'https://github.com/GeminiDRSoftware/GeminiCalMgr.git@release/1.1.x'
OBSDB_URL = R'https://github.com/GeminiDRSoftware/GeminiObsDB.git@release/1.0.x'

DRAGONS_BRANCH = 'master'
DRAGONS_LOCATION = 'DRAGONS/'

PATH = os.path.abspath(os.path.dirname(__file__))

# Define the following paths to be set as environment variables
NEW_ENV_VARIABLES = {    
    "MAROONX_LEGACY_TEST": Path("home/martin/Projects/MaroonX/legacy/maroonx_base/data2"),
    "MAROONX_DRAGONS_TEST": Path(PATH),
}


def check_dragons_version(session: nox.Session):
    """Check if dragons is the expected version."""
    with session.chdir(DRAGONS_LOCATION):
        result = session.run('git', 'branch', silent=True, external=True)

        match_ = re.match(r'^.*\s+(\w+)\s*.*$', result)

        if not match_:
            message = 'No DRAGONS branch found.'
            raise ValueError(message)

        branch_name = match_.group(1)

        if branch_name != DRAGONS_BRANCH:
            session.warn(f'Unexpected git branch: {branch_name} (not {DRAGONS_BRANCH})')

        else:
            session.log(f'Found correct branch: {branch_name}')

        result = session.run('git', 'fetch', '--dry-run', silent=True, external=True)

        if result:
            session.warn(
                f'Your DRAGONS version is not up-to-date.\n'
                f'Please check the latest version at:\n'
                f'    {DRAGONS_URL}\n'
                f'And, if you would like to update, run:\n\n'
                f'    git fetch && git pull\n\n'
                f' We strongly encourage you do this regularly in case of '
                f' important updates.'
            )

        else:
            session.log('DRAGONS is up to date!')


def install_dragons(session: nox.Session, python: Path | None = None):
    """Install dragons into the given session.

    If python is not None, it assumes it is a path to the
    correct python binary to use.
    """
    dragons_path = Path(DRAGONS_LOCATION)

    if not dragons_path.exists():
        # Clone dragons locally
        session.run(
            'git',
            'clone',
            '-b',
            DRAGONS_BRANCH,
            DRAGONS_URL,
            str(dragons_path),
            external=True,
        )

    check_dragons_version(session)

    if python:
        session.run(
            str(python),
            '-m',
            'pip',
            'install',
            '-e',
            str(dragons_path),
            external=True,
        )

        session.run(
            str(python),
            '-m',
            'pip',
            'install',
            f'git+{CALMGR_URL}',
            f'git+{OBSDB_URL}',
            external=True,
        )

        return

    session.install('-e', str(dragons_path))
    session.install(f'git+{CALMGR_URL}', f'git+{OBSDB_URL}')


def get_dependencies(session: nox.Session, only: str = '') -> list[str]:
    """Get dependencies from `pyproject.toml` using Poetry.

    Args:
        session: The nox session
        only: Optional dependency group name. If not provided, returns main
             dependencies. Valid values include "main", "test", "dev", or any
             custom group defined in pyproject.toml.

    Returns
    -------
        List of dependencies in the format "name==version"
    """
    lockfile_path = Path('poetry.lock')
    if not lockfile_path.exists():
        session.run('poetry', 'lock', external=True)

    only = only if only else 'main'

    cmd = ['poetry', 'show', '--top-level', '--only', only]

    result = session.run(
        *cmd,
        silent=True,
        external=True,
    )
    return _parse_dependencies(result)


def _parse_dependencies(result: str) -> list[str]:
    """Parse the output of poetry show to extract dependencies."""
    result = result.replace('(!)', ' ')
    dependencies = []
    for line in result.splitlines():
        if match := re.match(r'^\s*(\S+)\s+([\.0-9vV]+)\s*.*$', line):
            name = match.group(1)
            version = match.group(2)
            dependencies.append(f'{name}=={version}')
    return dependencies


# Development envs
@nox.session(venv_backend=None, python='3.12')
def devenv(session: nox.Session):
    """Create a development environment.

    This will perform the following steps:

    + Create a new virtual environment at ``venv/``
    + Install DRAGONS:
        + If DRAGONS does not exist locally, clone it.
        + Otherwise, perform a ``git fetch && git pull``
    + Install any other dependencies needed.
    """
    session.install('poetry', 'poetry-plugin-export')
    dependencies = get_dependencies(session, only='main,dev,test')

    env_name = 'mx_dev'
    session.run(
        'python3.12',
        '-m',
        'venv',
        'venv/',
        '--clear',
        '--upgrade-deps',
        '--prompt',
        env_name,
        external=True,
    )

    venv_loc = Path('venv').resolve()
    venv_python = venv_loc / 'bin' / 'python'

    # Install DRAGONS
    install_dragons(session, python=venv_python)

    # Install dependencies
    session.run(
        str(venv_python),
        '-m',
        'pip',
        'install',
        *dependencies,
        external=True,
    )

    # Install maroonxdr and maroonx_instruments
    session.run(
        str(venv_python),
        '-m',
        'pip',
        'install',
        '-e',
        '.',
        external=True,
    )

    # Add environment variables to the activate script
    venv_activate = venv_loc / 'bin' / 'activate'
        
    # Append environment variables to the activate script
    with open(venv_activate, 'a') as f:
        f.write('\n# Custom environment variables for MAROONXDR\n')
        for var_name, var_value in NEW_ENV_VARIABLES.items():
            f.write(f'export {var_name}="{str(var_value)}"\n')

    session.log(
        f'Successfully created virtual environment at {venv_loc}! '
        f'To activate your environment, run: \n'
        f'     source {venv_activate}\n'
        f'The following env variables are available: \n'
        f'     {list(NEW_ENV_VARIABLES.keys())}'
    )
    # session.notify('initialize_commit_hooks')


@nox.session(venv_backend=None, python='3.12')
def devconda(session: nox.Session):
    """Create a conda development environment."""
    session.install('poetry', 'poetry-plugin-export')
    dependencies = get_dependencies(session, only='main,dev,test')

    env_name = 'mx_devconda'
    session.run(
        'conda',
        'create',
        '--yes',
        '--force',
        '-n',
        env_name,
        '-c', 'http://astroconda.gemini.edu/public',
        '-c', 'conda-forge', 
        '-c', 'defaults',
        'python=3.12',
        external=True,
    )

    result = session.run('conda', 'info', '-e', silent=True, external=True)

    EXPECTED_COLUMNS = 2

    env_path = None

    for line in result.splitlines():
        info = line.split('#')[0]
        columns = info.split()

        if len(columns) != EXPECTED_COLUMNS:
            continue

        name, path = columns
        if name == env_name:
            env_path = Path(path)
            break

    if env_path is None:
        message = f'Could not find environment {env_name}'
        raise OSError(message)

    env_python = env_path / 'bin' / 'python'

    # Install dependencies from pyproject.toml without version constraints
    for dependency in dependencies:
        dep = dependency  # .split('==')[0]
        session.run(
            'conda', 'install', f'--name={env_name}', '--yes', dep, external=True
        )

    # Install DRAGONS conda dependencies as stated in its README
    session.run(
        'conda',
        'install',
        f'--name={env_name}',
        '--yes',
        '--no-update-deps',
        # '-c',
        # 'conda-forge',
        'astropy>=6',
        'astroquery',
        'matplotlib',
        'numpy<2',
        'psutil',
        'python-dateutil',
        'requests',
        'scikit-image',
        'scipy',
        'sextractor',
        'sqlalchemy>=2.0.0',
        'ds9',
        'gwcs>=0.15,<=0.22.1',
        'specutils',
        'sphinx',
        'sphinx_rtd_theme',
        'bokeh>=3',
        'holoviews',
        'cython',
        'future',
        'astroscrappy>=1.1',
        'fitsverify',
        'jsonschema',
        'imexam',
        external=True,
    )

    # Install DRAGONS
    install_dragons(session, python=env_python)

    # Install maroonxdr and maroonx_instruments
    session.run(
        str(env_python),
        '-m',
        'pip',
        'install',
        '-e',
        '.',
        external=True,
    )

    session.log('Conda environment generated, to activate run:')
    session.log(f'   conda activate {env_name}')


# Lint sessions
@nox.session
def lint(session: nox.Session):
    """Run linters."""
    session.install('ruff')
    session.run('ruff', 'check', '--fix')


@nox.session
def ruff_format(session: nox.Session):
    """Run formatters."""
    session.install('ruff')
    session.run('ruff', 'format')


@nox.session(python=False)
def initialize_commit_hooks(session: nox.Session):
    """Run pre-commit to install various hooks.

    The hooks are in `.pre-commit-config.yaml`.
    """
    # If not in a git repo, then pre-commit will fail. This should *not* be
    # considered an error, since the git repo should be set up by the person
    # running the script.
    if not Path('.git').exists():
        session.log(
            'Not in a git repository. Skipping pre-commit installation.'
            'If you meant to install pre-commit, please initialize '
            'this repository (`git init`).'
        )

        return

    session.install('pre-commit')

    # May be external --- doesn't need to be installed locally.
    session.run(
        'pre-commit',
        'install',
        '--install-hooks',
        '--hook-type=pre-commit',
        '--hook-type=commit-msg',
        external=True,
    )



@nox.session(python='3.12')
def unit_tests(session: nox.Session):
    """Run unit tests."""
    session.install('poetry', 'poetry-plugin-export')

    # Set environment variables that tests might need
    for var_name, var_value in NEW_ENV_VARIABLES.items():
        session.env[var_name] = str(var_value)

    # Install DRAGONS first
    install_dragons(session)

    # Install dependencies
    dependencies = get_dependencies(session, only='main,test')
    session.install(*dependencies)

    # Install maroonxdr and maroonx_instruments in editable mode
    session.install('-e', '.')

    # Run the tests with corrected paths
    test_args = [
        'pytest',
        'maroonxdr/maroonx/tests/',
        '--ignore=maroonxdr/maroonx/tests/complete',
        '--ignore=maroonxdr/maroonx/tests/regression',
        '-v',
        '--tb=no',
        '--rootdir=.',
    ]
    
    # Add any additional arguments passed via command line
    test_args.extend(session.posargs)
    session.run(*test_args)

@nox.session(python='3.12')
def regression_tests(session: nox.Session):
    """Run unit tests."""
    session.install('poetry', 'poetry-plugin-export')

    # Set environment variables that tests might need
    for var_name, var_value in NEW_ENV_VARIABLES.items():
        session.env[var_name] = str(var_value)

    # Install DRAGONS first
    install_dragons(session)

    # Install dependencies
    dependencies = get_dependencies(session, only='main,test')
    session.install(*dependencies)

    # Install maroonxdr and maroonx_instruments in editable mode
    session.install('-e', '.')

    # Run the tests with corrected paths
    test_args = [
        'pytest',
        'maroonxdr/maroonx/tests/regression',
        '-v',
        '--tb=no',
        '--rootdir=.',
    ]
    
    # Add any additional arguments passed via command line
    test_args.extend(session.posargs)
    session.run(*test_args)

@nox.session(python='3.12')
def complete_tests(session: nox.Session):
    """Run unit tests."""
    session.install('poetry', 'poetry-plugin-export')

    # Set environment variables that tests might need
    for var_name, var_value in NEW_ENV_VARIABLES.items():
        session.env[var_name] = str(var_value)

    # Install DRAGONS first
    install_dragons(session)

    # Install dependencies
    dependencies = get_dependencies(session, only='main,test')
    session.install(*dependencies)

    # Install maroonxdr and maroonx_instruments in editable mode
    session.install('-e', '.')

    # Run the tests with corrected paths
    completion_tests = [
        'maroonxdr/maroonx/tests/complete/bundle_test.py',
        'maroonxdr/maroonx/tests/complete/dark_test.py',
        'maroonxdr/maroonx/tests/complete/flat_test.py',
        'maroonxdr/maroonx/tests/complete/reduced_1D_test.py',
        'maroonxdr/maroonx/tests/complete/science_test.py',
    ]

    # Run each script using the session's Python
    for test_script in completion_tests:
        session.run('python', test_script) 

@nox.session(python='3.10')
def integration_tests(session: nox.Session):
    """Run integration tests."""
    message = f'{session.name} not configured.'
    raise NotImplementedError(message)


@nox.session(python='3.10')
def coverage(session: nox.Session):
    """Run tests with coverage reporting."""
    session.install('poetry', 'poetry-plugin-export')

    # Install DRAGONS
    install_dragons(session)

    # Install dependencies
    dependencies = get_dependencies(session, only='main,test')
    session.install(*dependencies)

    # Install maroonxdr and maroonx_instruments
    session.install('-e', '.')

    session.run('coverage', 'erase')
    session.run(
        'pytest',
        '-q',
        'maroonxdr/maroonx/tests/image/test_split_bundle.py',
        '--cov=maroonxdr/maroonx/',
        '--cov-append',
    )
    session.run('coverage', 'report', '--fail-under=80', '-m')


# Documentation
@nox.session(venv_backend='virtualenv')
def docs(session: nox.Session):
    """Build documentation using Sphinx."""
    session.install('poetry', 'poetry-plugin-export')
    
    # Install DRAGONS
    install_dragons(session)
    
    # Install dependencies
    dependencies = get_dependencies(session, only='main,docs')
    session.install(*dependencies)
    
    # Install maroonxdr and maroonx_instruments
    session.run(
        'pip',
        'install',
        '-e',
        '.',
        external=True,
    )
    
    # Create docs directories if they don't exist
    docs_source = Path('docs/source')
    docs_build = Path('docs/build/html')
    
    if not docs_source.exists():
        session.error(
            f"Documentation source directory {docs_source} does not exist. "
            "Please create it and add your Sphinx configuration."
        )
    
    # Build the documentation
    session.run(
        'sphinx-build', 
        '-b', 'html',
        '-W',  # Treat warnings as errors
        '--keep-going',  # Continue on errors to see all issues
        str(docs_source), 
        str(docs_build)
    )
    
    session.log(f"Documentation built successfully in {docs_build}")
    session.log(f"Open {docs_build}/index.html in your browser to view the docs")

@nox.session(venv_backend='virtualenv')
def docstyle(session: nox.Session):
    """Check docstring style using pydocstyle."""
    session.install('pydocstyle')
    session.run(
        'pydocstyle', 
        'maroonxdr/',
        '--convention=numpy',
        '--add-ignore=D100,D104'  # Ignore missing docstrings in modules and __init__.py
    )

# Building
@nox.session(venv_backend='virtualenv')
def build(session: nox.Session):
    """Build the project."""
    session.install('poetry')
    session.run('poetry', 'build', '--output=dist')
