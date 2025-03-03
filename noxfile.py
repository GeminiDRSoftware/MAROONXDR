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
nox.options.default_venv_backend = 'conda'
nox.options.error_on_external_run = True

# Dragons installation resources
DRAGONS_URL = R'https://github.com/GeminiDRSoftware/DRAGONS'
CALMGR_URL = R'https://github.com/GeminiDRSoftware/GeminiCalMgr.git@release/1.1.x'
OBSDB_URL = R'https://github.com/GeminiDRSoftware/GeminiObsDB.git@release/1.0.x'

DRAGONS_BRANCH = 'master'
DRAGONS_LOCATION = 'DRAGONS/'


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


@nox.session(venv_backend=None)
def devenv(session: nox.Session):
    """Create a development environment.

    This will perform the following steps:

    + Create a new virtual environment at ``venv/``
    + Install DRAGONS:
        + If DRAGONS does not exist locally, clone it.
        + Otherwise, perform a ``git fetch && git pull``
    + Install any other dependencies needed.
    """
    env_name = 'mx_dev'
    session.run(
        'python3.10',
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

    requirements_file = Path('requirements.txt')

    session.run(
        str(venv_python),
        '-m',
        'pip',
        'install',
        '-r',
        str(requirements_file),
        external=True,
    )

    venv_activate = venv_loc / 'bin' / 'activate'

    session.log(
        f'Successfully created virtual environment at {venv_loc}! '
        f'To activate your environment, run: \n'
        f'     source {venv_activate}\n'
    )

    session.notify('initialize_commit_hooks')


@nox.session(venv_backend=None)
def devconda(session: nox.Session):
    """Create a conda development environment."""
    env_name = 'mx_devconda'
    session.run(
        'conda',
        'create',
        '--yes',
        '--force',
        '-n',
        env_name,
        '-c',
        'conda-forge',
        'python=3.10',
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

    # Install DRAGONS
    install_dragons(session, python=env_python)

    session.log('Conda environment generated, to activate run:')
    session.log(f'   conda activate {env_name}')


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


@nox.session
def observer(session: nox.Session):
    """Run non-destructive checks."""
    session.install('ruff', 'pytest')

    # Check linting issues
    session.run('ruff', '-s', 'check')

    # Check tests
    session.run('pytest', 'maroonxdr/tests/', '--tb=no')


@nox.session(venv_backend=None)
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


# Testing
@nox.session(python='3.10')
def unit_tests(session: nox.Session):
    """Run unit tests."""
    session.install('pytest')
    session.run(
        'pytest', 'maroonxdr/tests/', 'maroonx_instruments/tests/', *session.posargs
    )


@nox.session(python='3.10')
def integration_tests(session: nox.Session):
    """Run integration tests."""
    message = f'{session.name} not configured.'
    raise NotImplementedError(message)


@nox.session
def coverage(session: nox.Session):
    """Run tests with coverage reporting."""
    session.install('coverage', 'pytest', 'pytest-cov')

    session.run('coverage', 'erase')
    session.run('pytest', '-q', 'maroonxdr/tests/', '--cov=maroonxdr/', '--cov-append')
    session.run('coverage', 'report', '--fail-under=80', '-m')


# Documentation
@nox.session(venv_backend='virtualenv')
def docs(session: nox.Session):
    """Build documentation using Sphinx."""
    session.install('sphinx', 'sphinx-rtd-theme', 'myst-parser')
    session.run('sphinx-build', '-b', 'html', 'docs/source', 'docs/build/html')


# Building
@nox.session(venv_backend='virtualenv')
def build(session: nox.Session):
    """Build the project."""
    session.install('poetry')
    session.run('poetry', 'build', '--output=dist')


# Code complexity
@nox.session(venv_backend='virtualenv')
def complexity(session: nox.Session):
    """Check code complexity metrics."""
    session.install('radon')
    session.run('radon', 'cc', 'maroonxdr/', '-a', '-n', 'C')
    session.run(
        'xenon',
        '--max-absolute',
        'C',
        '--max-modules',
        'C',
        '--max-average',
        'A',
        'maroonxdr/',
    )
