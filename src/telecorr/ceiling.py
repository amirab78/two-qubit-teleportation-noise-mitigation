"""The physical ceiling of a local unitary correction, and the correction that reaches it.

Roadmap: step 5 (incoherent noise and the ceiling).

Under the twirl (``telecorr.twirl``) the corrected leg keeps p_I = (1 + tr(R T)) / 4, so the best rotation R
maximises tr(R T) over SO(3): the orthogonal Procrustes problem, solved in closed form by the SVD of T with a
determinant fix (Kabsch). The maximum is the ceiling: no single-qubit unitary correction can do better.
"""
from __future__ import annotations

from typing import Mapping

import numpy as np

from .rotations import compose
from .twirl import f_avg_of


def procrustes(T: np.ndarray) -> np.ndarray:
    """Rotation R maximising tr(R T) (Kabsch / orthogonal Procrustes with det = +1)."""
    U, _, Wt = np.linalg.svd(T)
    d = np.sign(np.linalg.det(Wt.T @ U.T)) or 1.0
    return Wt.T @ np.diag([1.0, 1.0, d]) @ U.T


def procrustes_batch(T: np.ndarray) -> np.ndarray:
    """procrustes for a stack of matrices (n, 3, 3)."""
    U, _, Wt = np.linalg.svd(T)
    d = np.sign(np.linalg.det(np.swapaxes(Wt, 1, 2) @ np.swapaxes(U, 1, 2)))
    d[d == 0] = 1.0
    D = np.zeros_like(T); D[:, 0, 0] = 1; D[:, 1, 1] = 1; D[:, 2, 2] = d
    return np.swapaxes(Wt, 1, 2) @ D @ np.swapaxes(U, 1, 2)


def ceiling(T: Mapping[int, np.ndarray]) -> tuple[np.ndarray, float]:
    """(correction angles, F_avg) of the best message-independent local correction (exact, closed form)."""
    R = {q: procrustes(T[q]) for q in (2, 3)}
    return compose(None, R), f_avg_of(R, T)


def ceiling_pI(T: np.ndarray) -> float:
    """Largest p_I of one leg: (1 + sum of singular values with the det sign on the smallest) / 4."""
    s = np.linalg.svd(T, compute_uv=False)
    return float((1 + s[0] + s[1] + np.sign(np.linalg.det(T)) * s[2]) / 4)
