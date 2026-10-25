"""Counts-only process tomography of the channel legs: TOMO-9 and TOMO-12.

Roadmap: step 7 (efficient learning); run on real counts in step 8.

One probe round at setting S (a known rotation applied in the correction box, on both legs) teleports the three
probes |00>, |++>, |+i+i> and measures each in its own axis. Under the twirl the measured shrink factors are

    lambda[j, a] = sum_c S_j[a, c] T[c, a]      (linear in the 9 entries of T),

so 3 settings with 3 probes each (TOMO-9: S = I, R_cyc, R_cyc^2) fix T exactly, and the best correction is the
Procrustes rotation of the estimate. TOMO-12 uses the 4 settings {I, Rx(pi/2), Ry(pi/2), Rz(pi/2)}.
"""
from __future__ import annotations

from typing import Mapping, Sequence

import numpy as np

from .ceiling import procrustes
from .rotations import R_CYC, axis_rotation, compose, rotation_angle_axis, so3
from .twirl import f_avg

_R_SET = {None: np.eye(3), 0: so3(axis_rotation(0, np.pi / 2)), 1: so3(axis_rotation(1, np.pi / 2)),
          2: so3(axis_rotation(2, np.pi / 2))}
DESIGNS: dict[str, np.ndarray] = {"tomo12": np.stack([_R_SET[s] for s in (None, 0, 1, 2)]),
                                  "tomo9": np.stack([np.eye(3), R_CYC, R_CYC @ R_CYC])}
PROBE_AXES: tuple[str, str, str] = ("Z", "X", "Y")       # probe order within a setting: |00>, |++>, |+i+i>


def probe_labels(design: str = "tomo9") -> list[str]:
    """Binding labels of the probe circuits, in design order: ``"{design}_s{j}_{axis}"``."""
    return [f"{design}_s{j}_{axis}" for j in range(len(DESIGNS[design])) for axis in PROBE_AXES]


def ls_T(S: np.ndarray, lam: np.ndarray, ridge: float = 1e-9) -> np.ndarray:
    """Least-squares T of each leg from rounds at settings S (n, k, 3, 3) with results lam (n, k, 3)."""
    n, k = lam.shape[:2]
    A = np.zeros((n, k * 3, 9))
    for j in range(k):
        for a in range(3):
            A[:, 3 * j + a, [3 * c + a for c in range(3)]] = S[:, j, a, :]
    y = lam.reshape(n, -1)
    AtA = np.einsum("nri,nrj->nij", A, A) + ridge * np.eye(9)
    return np.linalg.solve(AtA, np.einsum("nri,nr->ni", A, y)[..., None])[..., 0].reshape(n, 3, 3)


def lambdas_from_counts(dists: Sequence[Mapping[str, float]], design: str = "tomo9") -> dict[int, np.ndarray]:
    """Shrink factors lambda[j, a] of each leg from Bob's frame-corrected counts of the probe circuits.
    Leg q2 is read on Bob's q4 (last key character), leg q3 on q5 (first character): lambda = 2 P(0) - 1."""
    k = len(DESIGNS[design])
    lam = {2: np.zeros((k, 3)), 3: np.zeros((k, 3))}
    for j in range(k):
        for i, axis in enumerate(PROBE_AXES):
            d = dists[3 * j + i]
            tot = sum(d.values())
            a = "XYZ".index(axis)
            lam[2][j, a] = 2 * sum(w for key, w in d.items() if key[-1] == "0") / tot - 1
            lam[3][j, a] = 2 * sum(w for key, w in d.items() if key[0] == "0") / tot - 1
    return lam


def estimate(dists: Sequence[Mapping[str, float]], design: str = "tomo9", x0: Sequence[float] | None = None) -> dict:
    """TOMO-9 / TOMO-12 from the probe counts (or exact probabilities), in the order of :func:`probe_labels`.

    Returns the T matrix of each leg (channel including the current correction ``x0``), the new correction angles
    (Procrustes rotation composed with ``x0``), the predicted F_avg before / after, the Pauli fingerprint, the
    coherent rotation found (angle, axis) and the singular values ("shrink") of T.
    """
    Rs = DESIGNS[design]
    lam = lambdas_from_counts(dists, design)
    T = {q: ls_T(Rs[None], lam[q][None], ridge=0.0)[0] for q in (2, 3)}
    Ropt = {q: procrustes(T[q]) for q in (2, 3)}
    x_new = compose(None if x0 is None else np.asarray(x0, float), Ropt)
    pI0 = {q: (1 + np.trace(T[q])) / 4 for q in (2, 3)}
    pI1 = {q: (1 + np.trace(Ropt[q] @ T[q])) / 4 for q in (2, 3)}
    out = {"design": design, "T": {q: T[q].tolist() for q in (2, 3)}, "x_new": [float(v) for v in x_new],
           "Favg_before": f_avg(pI0[2], pI0[3]), "Favg_after": f_avg(pI1[2], pI1[3]),
           "lambda": {q: lam[q].tolist() for q in (2, 3)}}
    for q in (2, 3):
        d = np.diag(T[q]); pI = (1 + d.sum()) / 4
        out[f"pauli_{q}"] = [float(pI)] + [float((1 + x) / 2 - pI) for x in d]
        ang, axis = rotation_angle_axis(Ropt[q].T)            # the rotation the correction undoes
        out[f"rotation_{q}"] = {"angle": ang, "axis": axis}
        out[f"shrink_{q}"] = np.linalg.svd(T[q], compute_uv=False).tolist()
    return out
