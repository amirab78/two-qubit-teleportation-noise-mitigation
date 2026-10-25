"""Crosstalk detector: are Bob's two output bits independent in the TOMO-9 probe circuits?

Roadmap: step 6 (correlated noise); run on real counts in step 8.

A local (product) channel keeps the two legs independent, so in every probe circuit the 2 x 2 table of Bob's
bits factorises. The detector is the Pearson chi-square test of independence summed over the 9 tables; a table
whose expected count falls below ``min_expected`` is skipped (the standard rule), each kept table adds one degree
of freedom. (A first version with a plug-in delta-method variance was anti-conservative and was replaced.)
"""
from __future__ import annotations

from typing import Mapping, Sequence

import numpy as np
from scipy.stats import chi2

from .paulis import KEYS4


def pearson_statistic(counts: Sequence[Mapping[str, float]], N: int, min_expected: float = 5) -> tuple[float, int]:
    """(X^2, degrees of freedom) over the given 2 x 2 tables ``{"b5b4": count}``; N = shots per circuit."""
    x2, dof = 0.0, 0
    for c in counts:
        n = np.array([c.get(k, 0) for k in KEYS4], float).reshape(2, 2)   # [b5][b4]
        E = np.outer(n.sum(1), n.sum(0)) / N
        if E.min() < min_expected:
            continue
        x2 += float(((n - E) ** 2 / E).sum())
        dof += 1
    return x2, dof


def detect(counts: Sequence[Mapping[str, float]], N: int, min_expected: float = 5) -> dict:
    """Detector result {"X2", "dof", "p"} (p = chi-square survival function, None without any usable table)."""
    x2, dof = pearson_statistic(counts, N, min_expected)
    return {"X2": float(x2), "dof": int(dof), "p": float(chi2.sf(x2, dof)) if dof else None}
