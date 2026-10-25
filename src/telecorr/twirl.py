"""The Pauli-twirl model of teleportation through a noisy channel.

Roadmap: step 2 (channel noise); the same model predicts the hardware in step 8.

Alice's four uniformly random outcome bits conjugate whatever happened on each channel leg by a random Pauli
(Pauli frame). Averaged over her outcomes, only the *diagonal* of each leg's Pauli transfer matrix reaches Bob:
with correction rotation R and channel T (3 x 3, unital part), Bob's leg is the Pauli channel

    lambda = diag(R T),   p_I = (1 + sum(lambda)) / 4,   p_a = (1 + lambda_a) / 2 - p_I.

For the two-qubit message the two legs act independently (q0 -> q4 through q2, q1 -> q5 through q3), so the
Haar-average fidelity of the teleported two-qubit state is F_avg = (4 p_I^(2) p_I^(3) + 1) / 5.
"""
from __future__ import annotations

from typing import Mapping

import numpy as np

from .paulis import PAULI


def leg_lambdas(R: np.ndarray, T: np.ndarray) -> np.ndarray:
    """Bloch shrink factors of Bob's Pauli channel on one leg: diag(R T)."""
    return np.einsum("ak,ka->a", R, T)


def twirled_map(R: np.ndarray, T: np.ndarray) -> np.ndarray:
    """The Bloch map that reaches Bob: diag(diag(R T))."""
    return np.diag(leg_lambdas(R, T))


def pauli_probs(lam: np.ndarray) -> np.ndarray:
    """(p_I, p_X, p_Y, p_Z) of a Pauli channel from its shrink factors."""
    pI = (1 + lam.sum()) / 4
    return np.r_[pI, (1 + lam) / 2 - pI]


def f_avg(pI2: float, pI3: float) -> float:
    """Haar-average fidelity of the teleported two-qubit state: (4 p_I2 p_I3 + 1) / 5."""
    return float((4 * pI2 * pI3 + 1) / 5)


def f_avg_of(R: Mapping[int, np.ndarray], T: Mapping[int, np.ndarray]) -> float:
    """F_avg for corrections R[q] and channels T[q] on the two legs (q = 2, 3)."""
    p = {q: (1 + np.trace(R[q] @ T[q])) / 4 for q in (2, 3)}
    return f_avg(p[2], p[3])


def leg_state_fidelity(psi: np.ndarray, M2: np.ndarray, M3: np.ndarray) -> float:
    """Bob's fidelity with the target ``psi`` when the message's two legs go through the unital Bloch maps
    M2 (q0 -> q4) and M3 (q1 -> q5). For teleportation pass the twirled maps diag(diag(R T))."""
    rho = np.outer(psi, psi.conj())
    P = [PAULI[a] for a in "IXYZ"]
    A = np.eye(4); A[1:, 1:] = np.asarray(M2)
    B = np.eye(4); B[1:, 1:] = np.asarray(M3)
    c = np.array([[np.real(np.trace(rho @ np.kron(P[b], P[a]))) for b in range(4)] for a in range(4)])  # a: q4, b: q5
    c2 = A @ c @ B.T
    out = sum(c2[a, b] * np.kron(P[b], P[a]) for a in range(4) for b in range(4)) / 4
    return float(np.real(psi.conj() @ out @ psi))


def message_fidelity(psi: np.ndarray, R: Mapping[int, np.ndarray], T: Mapping[int, np.ndarray]) -> float:
    """Exact teleportation fidelity of a pure two-qubit message with corrections R and channels T on the two legs."""
    return leg_state_fidelity(psi, twirled_map(R[2], T[2]), twirled_map(R[3], T[3]))
