"""Shared test fixtures / path setup.

Makes the package importable and regenerates the synthetic dataset into a
temporary directory so tests never depend on committed binary-ish artefacts.
Works under both `pytest` and `python -m unittest`.
"""

import os
import sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from examples.make_test_data import generate  # noqa: E402


def make_dataset(out_dir):
    """Generate the synthetic dataset into out_dir and return its info dict."""
    return generate(out_dir)
