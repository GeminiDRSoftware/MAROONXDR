"""Noxfile for the MAROONXDR package.

For a list of sessions and their descriptions, run:

.. code-block::
    nox -l

Please note that, to prevent undesirable execution, there are no default
sessions, so running ``nox`` in isolation will do nothing.
"""

import re
from pathlib import Path

import nox

nox.options.sessions = []
nox.options.error_on_external_run = True

# Dragons installation resources
DRAGONS_URL = R'https://github.com/GeminiDRSoftware/DRAGONS'
PYTEST_DRAGONS_URL = (
    R'git+https://github.com/GeminiDRSoftware/pytest_dragons.git@v1.0.0'
)
CALMGR_URL = (
    R'https://github.com/GeminiDRSoftware/GeminiCalMgr.git@release/1.1.x'  # deprecated
)
OBSDB_URL = (
    R'https://github.com/GeminiDRSoftware/GeminiObsDB.git@release/1.0.x'  # deprecated
)

# Meeting with Paul H. indicated that fits_storage is needed
# DRAGONS master (commit 22f4d9ff9) requires FitsStorage >= 3.4.0b1
FITSS_URL = R'https://github.com/GeminiDRSoftware/FitsStorage.git@v3.4.0b1'


DRAGONS_BRANCH = 'master'
DRAGONS_LOCATION = 'DRAGONS/'

PATH = Path(__file__).parent.resolve()

# Define the following paths to be set as environment variables
NEW_ENV_VARIABLES = {
    'DRAGONS_TEST': Path(PATH).parent / 'mx_test',
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
            f'git+{FITSS_URL}',
            external=True,
        )

        return

    session.install('-e', str(dragons_path))
    session.install(f'git+{FITSS_URL}')


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
    dependencies = get_dependencies(session, only='main,dev,docs,test')

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

    # Install pytest_dragons
    session.run(
        str(venv_python),
        '-m',
        'pip',
        'install',
        PYTEST_DRAGONS_URL,
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
    with venv_activate.open('a') as f:
        f.write('\n# Custom environment variables for MAROONXDR\n')
        for var_name, var_value in NEW_ENV_VARIABLES.items():
            f.write(f'export {var_name}="{var_value!s}"\n')

    session.log(
        f'Successfully created virtual environment at {venv_loc}! '
        f'To activate your environment, run: \n'
        f'     source {venv_activate}\n'
        f'The following env variables are available: \n'
        f'     {list(NEW_ENV_VARIABLES.keys())}'
    )


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
        '-c',
        'http://astroconda.gemini.edu/public',
        '-c',
        'conda-forge',
        '-c',
        'defaults',
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

    # Packages that are only available via pip (not in conda)
    # or that we prefer to install via pip
    pip_only_packages = {'barycorrpy', 'tables'}

    # Separate conda and pip dependencies
    conda_deps = []
    pip_deps = []
    for dependency in dependencies:
        package_name = dependency.split('==')[0].lower()
        if package_name in pip_only_packages:
            pip_deps.append(dependency)
        else:
            conda_deps.append(dependency)

    # Install conda dependencies
    for dep in conda_deps:
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

    # Install pip-only dependencies that aren't available in conda
    if pip_deps:
        session.run(
            str(env_python),
            '-m',
            'pip',
            'install',
            *pip_deps,
            external=True,
        )

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
def create_inputs(session: nox.Session):
    """Download and create input files for unit tests.
    Run this session to populate the input files before runing unit_tests.
    """
    session.install('poetry', 'poetry-plugin-export')

    # Set environment variables that tests might need
    for var_name, var_value in NEW_ENV_VARIABLES.items():
        session.env[var_name] = str(var_value)

    # Install DRAGONS first
    install_dragons(session)

    # Install dependencies
    dependencies = get_dependencies(session, only='main,test')
    session.install(*dependencies)

    # Install pytest_dragons
    session.install(f'{PYTEST_DRAGONS_URL}')

    # Install maroonxdr and maroonx_instruments in editable mode
    session.install('-e', '.')

    # Test modules that define create_inputs()
    create_inputs_scripts = [
        'maroonxdr/maroonx/tests/bundle/test_bundle.py',
        'maroonxdr/maroonx/tests/bundle/test_bundle_export.py',
        'maroonxdr/maroonx/tests/image/test_file_sorting.py',
        'maroonxdr/maroonx/tests/image/test_image_orientation_corrector.py',
        'maroonxdr/maroonx/tests/image/test_ND_filter_check.py',
        'maroonxdr/maroonx/tests/image/test_var.py',
    ]

    for script in create_inputs_scripts:
        session.run('python', script, '--create-inputs')


@nox.session(python='3.12')
def complete_tests(session: nox.Session):
    """Run complete end-to-end reduction tests."""
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

    # Install pytest_dragons
    session.install(f'{PYTEST_DRAGONS_URL}')

    # Run the tests with corrected paths
    completion_tests = [
        'maroonxdr/maroonx/tests/complete/bundle.py',
        'maroonxdr/maroonx/tests/complete/dark.py',
        'maroonxdr/maroonx/tests/complete/flat.py',
        'maroonxdr/maroonx/tests/complete/wavecal.py',
        'maroonxdr/maroonx/tests/complete/science.py',
    ]

    # Run each script using the session's Python
    for test_script in completion_tests:
        session.run('python', test_script, '--populate-inputs')


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

    # Install pytest_dragons
    session.install(f'{PYTEST_DRAGONS_URL}')

    # Install maroonxdr and maroonx_instruments in editable mode
    session.install('-e', '.')

    # Run the tests with corrected paths
    test_args = [
        'pytest',
        'maroonxdr/maroonx/tests/',
        '--ignore=maroonxdr/maroonx/tests/legacy_regression',
        '-v',
        '-rs',
        '--tb=short',
        '--rootdir=.',
    ]

    # Add any additional arguments passed via command line
    test_args.extend(session.posargs)
    session.run(*test_args)


@nox.session(python='3.12')
def legacy_regression_tests(session: nox.Session):
    """Run legacy pipeline regression tests."""
    session.install('poetry', 'poetry-plugin-export')

    # Set environment variables that tests might need
    for var_name, var_value in NEW_ENV_VARIABLES.items():
        session.env[var_name] = str(var_value)

    # Install DRAGONS first
    install_dragons(session)

    # Install dependencies
    dependencies = get_dependencies(session, only='main,test')
    session.install(*dependencies)

    # Install pytest_dragons
    session.install(f'{PYTEST_DRAGONS_URL}')

    # Install maroonxdr and maroonx_instruments in editable mode
    session.install('-e', '.')

    # Run the tests with corrected paths
    test_args = [
        'pytest',
        'maroonxdr/maroonx/tests/legacy_regression',
        '-v',
        '--tb=no',
        '--rootdir=.',
    ]

    # Add any additional arguments passed via command line
    test_args.extend(session.posargs)
    session.run(*test_args)


@nox.session(python='3.12')
def regression_tests(session: nox.Session):
    """Run DRAGONS-style regression tests with inputs/refs comparison."""
    session.install('poetry', 'poetry-plugin-export')

    # Set environment variables that tests might need
    for var_name, var_value in NEW_ENV_VARIABLES.items():
        session.env[var_name] = str(var_value)

    # Install DRAGONS first
    install_dragons(session)

    # Install dependencies
    dependencies = get_dependencies(session, only='main,test')
    session.install(*dependencies)

    # Install pytest_dragons
    session.install(f'{PYTEST_DRAGONS_URL}')

    # Install maroonxdr and maroonx_instruments in editable mode
    session.install('-e', '.')

    # Run the tests with corrected paths
    test_args = [
        'pytest',
        'maroonxdr/maroonx/tests/regression',
        '-v',
        '--tb=short',
        '--rootdir=.',
    ]

    # Add any additional arguments passed via command line
    test_args.extend(session.posargs)
    session.run(*test_args)


@nox.session(python='3.12')
def integration_tests(session: nox.Session):
    """Run integration tests."""
    message = f'{session.name} not configured.'
    raise NotImplementedError(message)


@nox.session(python='3.12')
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

    # Install pytest_dragons
    session.install(f'{PYTEST_DRAGONS_URL}')

    session.run('coverage', 'erase')
    session.run(
        'pytest',
        '-q',
        'maroonxdr/maroonx/tests/',
        '--ignore=maroonxdr/maroonx/tests/legacy_regression',
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

    # Define documentation directories
    doc_dir = Path('doc')
    user_source = doc_dir / 'usermanuals' / 'MAROONXDR_UserManual'
    user_build = user_source / 'build' / 'html'
    prog_source = doc_dir / 'progmanuals' / 'MAROONXDR_ProgManual'
    prog_build = prog_source / 'build' / 'html'
    tutorial_source = doc_dir / 'tutorials' / 'MAROONXDR_Tutorial'
    tutorial_build = tutorial_source / 'build' / 'html'

    # Build user manual
    session.log('Building user manual...')
    session.run(
        'sphinx-build', '-M', 'html', str(user_source), str(user_source / 'build')
    )

    # Build programmer manual
    session.log('Building programmer manual...')
    session.run(
        'sphinx-build', '-M', 'html', str(prog_source), str(prog_source / 'build')
    )

    # Build tutorial
    session.log('Building tutorial...')
    session.run(
        'sphinx-build',
        '-M',
        'html',
        str(tutorial_source),
        str(tutorial_source / 'build'),
    )

    session.log('Documentation built successfully!')
    session.log(f'User manual: {user_build}/index.html')
    session.log(f'Programmer manual: {prog_build}/index.html')
    session.log(f'Tutorial: {tutorial_build}/index.html')


@nox.session(venv_backend='virtualenv')
def docstyle(session: nox.Session):
    """Check docstring style using pydocstyle."""
    session.install('pydocstyle')
    session.run(
        'pydocstyle',
        'maroonxdr/',
        '--convention=numpy',
        '--add-ignore=D100,D104',  # Ignore missing docstrings in modules and __init__.py
    )
