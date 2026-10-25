"""Saved results of the research project that the notebooks and tests compare with (``data/reference``).

Utility: used by notebooks 1-7 and the regression tests.
"""
from __future__ import annotations

import json
import os
from functools import lru_cache

from .hardware.data import data_dir


def simulation(name: str) -> dict:
    """A saved simulation result, e.g. ``simulation("6a_nominal.json")``."""
    with open(os.path.join(data_dir(), "reference", "simulation", name), encoding="utf8") as fh:
        return json.load(fh)


def path(name: str) -> str:
    return os.path.join(data_dir(), "reference", "simulation", name)


@lru_cache(maxsize=1)
def _numbers() -> dict:
    with open(os.path.join(data_dir(), "reference", "final_numbers.json"), encoding="utf8") as fh:
        return json.load(fh)["numbers"]


def number(key: str):
    """A value of the project's final number table (``final_numbers.json``), e.g. ``number("n0003")``."""
    return _numbers()[key]["value"]


def printed(key: str) -> str:
    """The value of the final number table as the study printed it, e.g. ``printed("n0017")`` -> "9.3·10⁻⁶"."""
    return _numbers()[key]["printed"]
