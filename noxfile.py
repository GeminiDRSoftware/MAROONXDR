"""Noxfile for the MAROONXDR package.

For a list of sessions and their descriptions, run:

.. code-block::
    nox -l

Please note that, to prevent undesirable execution, there are no default
sessions, so running ``nox`` in isolation will do nothing.
"""

import datetime
import re
import shutil
import tomllib
import zipfile
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


GIT_ONLY_PACKAGES = {'pytest-dragons'}


def _parse_dependencies(result: str) -> list[str]:
    """Parse the output of poetry show to extract dependencies."""
    result = result.replace('(!)', ' ')
    dependencies = []
    for line in result.splitlines():
        if match := re.match(r'^\s*(\S+)\s+([\.0-9vV]+)\s*.*$', line):
            name = match.group(1)
            version = match.group(2)
            if name in GIT_ONLY_PACKAGES:
                continue
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
def preprocess(session: nox.Session):
    """Run the full blessed reduction chain into $DRAGONS_TEST/preprocessed_files/.

    Runs the five scripts under ``maroonxdr/maroonx/tests/preprocess/`` in
    order (bundle, dark, flat, wavecal, science). Their outputs are the
    blessed products that ``create_inputs`` nox session stages into the
    ``inputs/`` and ``refs/`` directories. This step is slow.
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

    # Install maroonxdr and maroonx_instruments in editable mode
    session.install('-e', '.')

    # Install pytest_dragons
    session.install(f'{PYTEST_DRAGONS_URL}')

    preprocess_scripts = [
        'maroonxdr/maroonx/tests/preprocess/bundle.py',
        'maroonxdr/maroonx/tests/preprocess/dark.py',
        'maroonxdr/maroonx/tests/preprocess/flat.py',
        'maroonxdr/maroonx/tests/preprocess/wavecal.py',
        'maroonxdr/maroonx/tests/preprocess/science.py',
    ]

    for script in preprocess_scripts:
        session.run('python', script, *session.posargs)


@nox.session(python='3.12')
def create_inputs(session: nox.Session):
    """Stage per-module test inputs and references.

    Run ``preprocess`` nox session before this.
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

    # Test modules that define staging hooks.
    create_inputs_modules = [
        'maroonx_instruments.maroonx.tests.test_calibration',
        'maroonxdr.maroonx.tests.image.test_remove_straylight',
        'maroonxdr.maroonx.tests.image.test_stack_flats',
        'maroonxdr.maroonx.tests.echelle_extraction.test_get_peaks_and_polynomials',
        'maroonxdr.maroonx.tests.echelle_extraction.test_fit_and_apply_etalon_wls',
        'maroonxdr.maroonx.tests.echelle_extraction.test_optimal_extraction',
        'maroonxdr.maroonx.tests.echelle_extraction.test_apply_wavelength_solution',
        'maroonxdr.maroonx.tests.echelle_extraction.test_box_extraction',
        'maroonxdr.maroonx.tests.echelle_extraction.test_combine_fibers',
        'maroonxdr.maroonx.tests.echelle_extraction.test_barycentric_correction',
    ]

    # Modules whose references are separate files in refs/ (the others embed
    # the reference in the input file itself).
    create_refs_modules = [
        'maroonxdr.maroonx.tests.echelle_extraction.test_optimal_extraction',
        'maroonxdr.maroonx.tests.echelle_extraction.test_box_extraction',
    ]

    for module in create_inputs_modules:
        session.run('python', '-m', module, '--create-inputs')

    for module in create_refs_modules:
        session.run('python', '-m', module, '--create-refs')


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
        'maroonx_instruments/maroonx/tests/',
        '-m', 'not preprocessed_data',
        '-v',
        '-rs',
        '--tb=short',
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
        'maroonxdr/maroonx/tests/',
        '-m', 'regression',
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

    # Set environment variables that tests might need
    for var_name, var_value in NEW_ENV_VARIABLES.items():
        session.env[var_name] = str(var_value)

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
        'maroonx_instruments/maroonx/tests/',
        '--cov=maroonxdr/maroonx/',
        '--cov=maroonx_instruments/maroonx/',
        '--cov-append',
    )
    session.run('coverage', 'report', '-m')


@nox.session(venv_backend='none')
def package_test_data(session: nox.Session):
    """Package the staged test data for sharing with other developers.

    Zips the per-module ``inputs/`` and ``refs/`` trees under
    ``$DRAGONS_TEST`` (the ``maroonxdr/`` and ``maroonx_instruments/``
    directories), excluding log files. The archive is written next to the
    data root. Recipients unzip it into their own ``$DRAGONS_TEST``
    directory and run the ``unit_tests`` and ``regression_tests``
    sessions; neither ``preprocess`` nor ``create_inputs`` is needed on
    their side.
    """
    dragons_test = Path(NEW_ENV_VARIABLES['DRAGONS_TEST'])
    data_trees = ['maroonxdr', 'maroonx_instruments']

    present = [d for d in data_trees if (dragons_test / d).is_dir()]
    if not present:
        message = f'No staged test data found under {dragons_test}'
        raise FileNotFoundError(message)
    for missing in sorted(set(data_trees) - set(present)):
        session.warn(f'{dragons_test / missing} not found; not packaged.')

    # build output filename
    stamp = datetime.datetime.now(tz=datetime.UTC).strftime('%Y%m%d')
    pyproject = tomllib.loads((PATH / 'pyproject.toml').read_text())
    version = pyproject['tool']['poetry']['version']
    archive = dragons_test.parent / f'mx_test_{version}_{stamp}.zip'

    n_files = 0
    with zipfile.ZipFile(archive, 'w', compression=zipfile.ZIP_DEFLATED) as zf:
        for tree in present:
            for path in sorted((dragons_test / tree).rglob('*')):
                if path.is_dir() or path.suffix == '.log':
                    continue
                zf.write(path, path.relative_to(dragons_test))
                n_files += 1

    size_gb = archive.stat().st_size / 1024**3
    session.log(f'Packaged {n_files} files into {archive} ({size_gb:.1f} GB)')


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
        # Ignore missing docstrings in modules and __init__.py
        '--add-ignore=D100,D104',
    )


# Primitives to exclude from the user-manual auto-generated reference:
# deprecated "Old" variants and the regression-only legacy patch.
_USERMANUAL_SKIP_PRIMITIVES = (
    'stackDarksOld',
    'stackFlatsOld',
    'removeStrayLight_legacyPatch',
)


def _regen_usermanual_autodoc(session: nox.Session, source_dir: Path):
    """Regenerate the user-manual primitive and recipe .rst fragments.

    Runs the ``generate_primdoc.py`` and ``generate_recipedoc.py`` utility
    scripts against the live ``maroonxdr`` code, writing fragments into the
    ``generated_doc/`` directories that the hand-written pages ``.. include``.
    The generators import the pipeline directly, so an ``mx_dev`` environment
    must be active (``venv_backend='none'``).

    The deprecated/legacy-patch primitive fragments listed in
    ``_USERMANUAL_SKIP_PRIMITIVES`` are pruned after generation so they do not
    appear in the reference.
    """
    scripts = Path('doc') / 'usermanuals' / 'utility_scripts'
    prim_gen = source_dir / 'primitives' / 'generated_doc'
    rec_gen = source_dir / 'recipes' / 'generated_doc'

    prim_gen.mkdir(parents=True, exist_ok=True)
    rec_gen.mkdir(parents=True, exist_ok=True)

    session.log('Regenerating primitive fragments...')
    session.run(
        'python', str(scripts / 'generate_primdoc.py'),
        'maroonx', '-d', str(prim_gen),
        external=True,
    )

    session.log('Regenerating recipe fragments (sq, qa)...')
    for context in ('sq', 'qa'):
        session.run(
            'python', str(scripts / 'generate_recipedoc.py'),
            'maroonx', '--context', context, '-d', str(rec_gen),
            external=True,
        )

    # Prune deprecated / legacy-patch primitive fragments.
    for primitive in _USERMANUAL_SKIP_PRIMITIVES:
        for fragment in prim_gen.glob(f'*.MAROONX.{primitive}_*.rst'):
            session.log(f'Pruning {fragment.name}')
            fragment.unlink()


def _build_manual(session: nox.Session, source_dir: Path, name: str):
    """Build a manual as HTML (default) or PDF.

    Pass ``--pdf`` in the nox posargs to build a PDF via LaTeX instead.
    """
    build_pdf = '--pdf' in session.posargs
    target = 'latexpdf' if build_pdf else 'html'
    build_dir = source_dir / 'build'

    if build_dir.exists():
        session.log(f'Cleaning {build_dir}')
        shutil.rmtree(build_dir)

    session.log(f'Building {name} ({target})...')
    session.run(
        'sphinx-build', '-M', target, str(source_dir), str(build_dir),
        external=True,
    )

    if build_pdf:
        pdfs = sorted((build_dir / 'latex').glob('*.pdf'))
        if pdfs:
            uri = pdfs[0].resolve().as_uri()
            session.log(f'{name} PDF: {uri}')
        else:
            session.log(
                f'{name} PDF build finished but no .pdf found under '
                f'{build_dir}/latex/ - check the build output above.'
            )
    else:
        uri = (build_dir / 'html' / 'index.html').resolve().as_uri()
        session.log(f'{name} HTML: {uri}')


@nox.session(venv_backend='none')
def usermanual(session: nox.Session):
    """Build the MaroonX User Manual.

    Defaults to HTML. Pass ``--pdf`` to build a PDF instead:
        nox -s usermanual -- --pdf

    Pass ``--regen`` to regenerate the primitive and recipe reference
    fragments from the live docstrings before building:
        nox -s usermanual -- --regen

    Requires ``mx_dev`` activated (Sphinx + DRAGONS on PATH).
    """
    source_dir = Path('doc/usermanuals/MAROONXDR_UserManual')

    if '--regen' in session.posargs:
        _regen_usermanual_autodoc(session, source_dir)

    _build_manual(session, source_dir, 'User Manual')


@nox.session(venv_backend='none')
def progmanual(session: nox.Session):
    """Build the MaroonX Programmer Manual.

    Defaults to HTML. Pass ``--pdf`` to build a PDF instead:
        nox -s progmanual -- --pdf

    Requires ``mx_dev`` activated (Sphinx + DRAGONS on PATH).
    """
    _build_manual(
        session,
        Path('doc/progmanuals/MAROONXDR_ProgManual'),
        'Programmer Manual',
    )


@nox.session(venv_backend='none')
def tutorial(session: nox.Session):
    """Build the MaroonX Tutorial.

    Defaults to HTML. Pass ``--pdf`` to build a PDF instead:
        nox -s tutorial -- --pdf

    Requires ``mx_dev`` activated (Sphinx + DRAGONS on PATH).
    """
    _build_manual(
        session,
        Path('doc/tutorials/MAROONXDR_Tutorial'),
        'Tutorial',
    )