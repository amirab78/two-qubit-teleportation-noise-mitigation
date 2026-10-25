"""Pauli matrices and small helpers shared by every module.

Roadmap: shared by every step.

Conventions (used throughout the package):

* qubit order is qiskit's little-endian order: a two-qubit operator on (q4, q5) is written ``kron(A_q5, A_q4)``;
* outcome keys of Bob's two bits are strings ``"b5b4"`` (qiskit order), e.g. ``"01"`` = q5 read 0, q4 read 1;
* a Bloch vector is ordered (x, y, z); a 3 x 3 matrix acting on Bloch vectors is a "T matrix" or rotation.
"""
from __future__ import annotations

import numpy as np

I2: np.ndarray = np.eye(2, dtype=complex)
X: np.ndarray = np.array([[0, 1], [1, 0]], dtype=complex)
Y: np.ndarray = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z: np.ndarray = np.array([[1, 0], [0, -1]], dtype=complex)

PAULI: dict[str, np.ndarray] = {"I": I2, "X": X, "Y": Y, "Z": Z}
SIGMA: tuple[np.ndarray, np.ndarray, np.ndarray] = (X, Y, Z)          # the Bloch-vector generators
KEYS4: tuple[str, str, str, str] = ("00", "01", "10", "11")          # Bob's outcome keys "b5b4"
TOMO_BASES: tuple[str, ...] = tuple(a + b for a in "XYZ" for b in "XYZ")  # (q4 basis, q5 basis)


def bloch(rho: np.ndarray) -> np.ndarray:
    """Bloch vector (x, y, z) of a one-qubit density matrix."""
    return np.array([np.real(np.trace(rho @ s)) for s in SIGMA])


def pauli_2q(p4: str, p5: str) -> np.ndarray:
    """The two-qubit Pauli product sigma_p5 (x) sigma_p4 (qiskit order)."""
    return np.kron(PAULI[p5], PAULI[p4])


def correlation_matrix(psi: np.ndarray) -> np.ndarray:
    """c[u4, u5] = <psi| sigma_u5 (x) sigma_u4 |psi> for u = I, X, Y, Z (index 0..3)."""
    L = "IXYZ"
    return np.array([[np.real(np.conj(psi) @ pauli_2q(L[a], L[b]) @ psi) for b in range(4)] for a in range(4)])
