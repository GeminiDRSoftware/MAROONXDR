"""Nox sessions for various automations."""

import re
from pathlib import Path

import nox

nox.options.sessions = ["lint", "ruff_format"]


def get_dependencies(session: nox.Session) -> list[str]:
    """Get dependencies from `pyproject.toml` using Poetry."""
    lockfile_path = Path("poetry.lock")

    if not lockfile_path.exists():
        session.run("poetry", "lock", external=True)

    result = session.run(
        "poetry",
        "show",
        "--top-level",
        silent=True,
        external=True,
    )

    result = result.replace("(!)", "   ")

    dependencies = []
    for line in result.splitlines():
        if match := re.match(r"^\s*(\S+)\s+([\.0-9vV]+)\s*.*$", line):
            name = match.group(1)
            version = match.group(2)
            dependencies.append(f"{name}=={version}")

    return dependencies


@nox.session
def create_venv(session: nox.Session):
    """Create a new virtual environment for development."""
    session.install("poetry", "poetry-plugin-export", "virtualenv")
    dependencies = get_dependencies(session)

    session.run(
        "virtualenv",
        "venv",
        "--prompt=(MAROONXDR)",
    )

    # Install dependencies
    session.run("venv/bin/pip", "install", *dependencies, external=True)
    session.log("Dependencies installed. Listing them below:")
    session.run("venv/bin/pip", "freeze", external=True)

    # Initialize pre-commit hooks
    session.notify("initialize_commit_hooks")


@nox.session(venv_backend="conda")
def create_conda_env(session: nox.Session):
    """Create a new conda environment for development."""
    session.install("poetry", "poetry-plugin-export")

    # Get the dependencies directly using get_dependencies
    dependencies = get_dependencies(session)

    # Create the Conda environment
    session.run(
        "conda",
        "create",
        "--force",
        "--name=mx_dragons",
        "--yes",
        "python=3.9",  # You can specify the version of Python if needed
    )

    # Install dependencies using conda for each dependency
    for dep in dependencies:
        session.run(
            "conda", "install", "--name=mx_dragons", "--yes", dep, external=True
        )

    # Initialize pre-commit hooks
    session.notify("initialize_commit_hooks")


@nox.session
def lint(session: nox.Session):
    """Run linters."""
    session.install("ruff")
    session.run("ruff", "check", "--fix")


@nox.session
def ruff_format(session: nox.Session):
    """Run formatters."""
    session.install("ruff")
    session.run("ruff", "format")


@nox.session
def initialize_commit_hooks(session: nox.Session):
    """Run pre-commit to install various hooks.

    The hooks are in `.pre-commit-config.yaml`.
    """
    # If not in a git repo, then pre-commit will fail. This should *not* be
    # considered an error, since the git repo should be set up by the person
    # running the script.
    if not Path(".git").exists():
        session.log(
            "Not in a git repository. Skipping pre-commit installation."
            "If you meant to install pre-commit, please initialize "
            "this repository (`git init`)."
        )

        return

    session.install("pre-commit")

    # May be external --- doesn't need to be installed locally.
    session.run(
        "pre-commit",
        "install",
        "--install-hooks",
        "--hook-type=pre-commit",
        "--hook-type=commit-msg",
        external=True,
    )


# Testing
@nox.session(python=["3.10", "3.11", "3.12"])
def unit_tests(session: nox.Session):
    """Run unit tests."""
    session.install("pytest")
    session.run("pytest", "tests/unit", *session.posargs)


@nox.session(python=["3.10", "3.11", "3.12"])
def integration_tests(session: nox.Session):
    """Run integration tests."""
    session.install("pytest")
    session.run("pytest", "tests/integration", *session.posargs)


@nox.session(python=["3.10", "3.11", "3.12"])
def build_tests(session: nox.Session):
    """Run build tests."""


# Building
@nox.session
def build(session: nox.Session):
    """Build the project."""
    session.install("poetry")
    session.run("poetry", "build", "--output=dist")
