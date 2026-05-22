"""Shared test utilities."""

import os
from contextlib import contextmanager
from pathlib import Path


@contextmanager
def change_cwd_context(target_dir):
    """
    Context manager to temporarily change working directory.

    Parameters
    ----------
    target_dir : str or Path
        Directory to change to

    Yields
    ------
    Path
        Path object of the target directory
    """
    original_dir = Path.cwd()
    target_path = Path(target_dir)
    try:
        os.chdir(target_path)
        yield target_path
    finally:
        os.chdir(original_dir)
