"""Single-qubit quantum channels: Kraus, Choi and Pauli-transfer forms, and the noise used in the study.

Roadmap: step 2 (coherent channel noise) and step 5 (incoherent noise).

Representations of a channel N(rho) = sum_k K_k rho K_k^dag on one qubit:

* Kraus operators ``[K_k]``;
* Choi matrix J = sum_ij |i><j| (x) N(|i><j|) (trace 2; positive <=> completely positive);
* Pauli transfer matrix (PTM) R_ab = tr(s_a N(s_b)) / 2 (a, b = I, X, Y, Z); its lower-right 3 x 3 block is the
  matrix **T** that acts on Bloch vectors and is all that teleportation lets through (``telecorr.twirl``).

The channel of each channel qubit in the study is "coherent rotation, then incoherent part": the rotation is
exp(-i angle/2 n.sigma); the incoherent part is amplitude damping (gamma), then dephasing (lambda), then
depolarizing (d). Three families scale with one strength p (q3 gets 1.5 p): D (depolarizing p), A (amplitude
damping p) and M (amplitude damping p/2 + dephasing p).
"""
from __future__ import annotations

from typing import Mapping, Sequence

import numpy as np

from .paulis import I2, SIGMA
from .rotations import NOMINAL_NOISE, rotation, so3

PAULIS4 = (I2,) + SIGMA
FAMILIES: dict[str, tuple[str, callable]] = {
    "D": ("depolarizing", lambda p: {"depol": p}),
    "A": ("amplitude damping", lambda p: {"amp": p}),
    "M": ("amplitude damping p/2 + dephasing p", lambda p: {"amp": p / 2, "phase": p}),
}
RATIO = {2: 1.0, 3: 1.5}                       # q3 gets 1.5 x the strength of q2
#: the "hard case" of Stage 4 (on which a single-start calibration got stuck), drawn as in Stage 4
HARD_CASE_SEED = 11


# ---------------------------------------------------------------- representations
def choi(kraus: Sequence[np.ndarray]) -> np.ndarray:
    """Choi matrix sum_ij |i><j| (x) N(|i><j|) (d^2 x d^2)."""
    d = kraus[0].shape[0]
    J = np.zeros((d * d, d * d), complex)
    for i in range(d):
        for j in range(d):
            E = np.zeros((d, d)); E[i, j] = 1
            J += np.kron(E, sum(K @ E @ K.conj().T for K in kraus))
    return J


def kraus_from_choi(J: np.ndarray, tol: float = 1e-12) -> list[np.ndarray]:
    """A minimal Kraus set from a Choi matrix (eigen-decomposition)."""
    d = int(round(np.sqrt(J.shape[0])))
    w, v = np.linalg.eigh((J + J.conj().T) / 2)
    return [np.sqrt(wk) * v[:, k].reshape(d, d).T for k, wk in enumerate(w) if wk > tol]


def ptm(kraus: Sequence[np.ndarray]) -> np.ndarray:
    """Full 4 x 4 Pauli transfer matrix R_ab = tr(s_a N(s_b)) / 2 of a one-qubit channel."""
    return np.array([[0.5 * np.real(sum(np.trace(Pa @ K @ Pb @ K.conj().T) for K in kraus)) for Pb in PAULIS4]
                     for Pa in PAULIS4])


def T_matrix(kraus: Sequence[np.ndarray]) -> np.ndarray:
    """The unital 3 x 3 block T of the PTM (what acts on Bloch vectors)."""
    return ptm(kraus)[1:, 1:]


def is_cptp(kraus: Sequence[np.ndarray], tol: float = 1e-12) -> bool:
    """Completely positive (Kraus form) and trace preserving: sum K^dag K = I."""
    d = kraus[0].shape[0]
    return bool(np.abs(sum(K.conj().T @ K for K in kraus) - np.eye(d)).max() < tol)


def is_unital(kraus: Sequence[np.ndarray], tol: float = 1e-12) -> bool:
    d = kraus[0].shape[0]
    return bool(np.abs(sum(K @ K.conj().T for K in kraus) - np.eye(d)).max() < tol)


# ---------------------------------------------------------------- the noise of the study
def incoherent_kraus(amp: float = 0.0, phase: float = 0.0, depol: float = 0.0) -> list[np.ndarray]:
    """Amplitude damping (gamma = amp), then dephasing (lambda = phase), then depolarizing
    (rho -> (1 - d) rho + d I/2), as a minimal Kraus set."""
    ad = [np.array([[1, 0], [0, np.sqrt(1 - amp)]]), np.array([[0, np.sqrt(amp)], [0, 0]])]
    pd = [np.array([[1, 0], [0, np.sqrt(1 - phase)]]), np.array([[0, 0], [0, np.sqrt(phase)]])]
    dp = [np.sqrt(1 - 3 * depol / 4) * I2] + [np.sqrt(depol / 4) * P for P in SIGMA]
    ks = [(D @ P @ A).astype(complex) for A in ad for P in pd for D in dp]
    return kraus_from_choi(choi(ks))


def incoherent_T(amp: float = 0.0, phase: float = 0.0, depol: float = 0.0) -> np.ndarray:
    """T of the incoherent part in closed form: diag(a, a, (1 - gamma)(1 - d)), a = sqrt((1-gamma)(1-lambda)) (1-d)."""
    a = np.sqrt((1 - amp) * (1 - phase)) * (1 - depol)
    return np.diag([a, a, (1 - amp) * (1 - depol)])


def leg_parameters(family: str | None, p: float, q: int) -> dict:
    """Incoherent parameters of channel qubit q (strength p on q2, 1.5 p on q3)."""
    cfg = {"amp": 0.0, "phase": 0.0, "depol": 0.0}
    if family is not None and p > 0:
        cfg.update(FAMILIES[family][1](RATIO[q] * p))
    return cfg


def coherent_unitaries(coherent: Mapping[int, Mapping] | None = None, strength: float = 1.0) -> dict[int, np.ndarray]:
    """Rotation of each channel qubit from {q: {"axis", "angle"}} (default: the Stage 2 noise N1), angles x strength."""
    coherent = NOMINAL_NOISE if coherent is None else coherent
    return {q: rotation(c["axis"], strength * c["angle"]) for q, c in coherent.items()}


def leg_kraus(coherent: Mapping[int, Mapping] | None = None, family: str | None = None, p: float = 0.0,
              strength: float = 1.0) -> dict[int, list[np.ndarray]]:
    """Kraus operators of each channel qubit: the coherent rotation first, then the incoherent part."""
    U = coherent_unitaries(coherent, strength)
    return {q: [K @ U[q] for K in incoherent_kraus(**leg_parameters(family, p, q))] for q in (2, 3)}


def leg_T(coherent: Mapping[int, Mapping] | None = None, family: str | None = None, p: float = 0.0,
          strength: float = 1.0) -> dict[int, np.ndarray]:
    """T of each channel qubit, T = D_incoherent R_rotation (closed form, = T_matrix(leg_kraus(...)))."""
    U = coherent_unitaries(coherent, strength)
    return {q: incoherent_T(**leg_parameters(family, p, q)) @ so3(U[q]) for q in (2, 3)}


def pair_kraus(legs: Mapping[int, Sequence[np.ndarray]]) -> list[np.ndarray]:
    """Two-qubit Kraus operators (order q3 (x) q2) of independent channels on q2 and q3."""
    return [np.kron(b, a) for b in legs[3] for a in legs[2]]


def random_noise(rng: np.random.Generator, max_angle: float = np.pi) -> dict[int, dict]:
    """Random coherent noise: on each channel qubit an axis uniform on the sphere and an angle in [0, max_angle]."""
    return {q: {"axis": tuple(rng.normal(size=3)), "angle": float(rng.uniform(0, max_angle))} for q in (2, 3)}


def hard_case() -> dict[int, dict]:
    """The Stage 4 "hard case": the 4th draw of rng(11) with angles in [0.3, 3.0]."""
    rng = np.random.default_rng(HARD_CASE_SEED)
    for _ in range(4):
        cfg = {q: {"axis": tuple(rng.normal(size=3)), "angle": float(rng.uniform(0.3, 3.0))} for q in (2, 3)}
    return cfg
