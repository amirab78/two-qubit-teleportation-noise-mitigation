"""Single-qubit unitaries, their Bloch rotations and the correction angles.

Roadmap: step 2 (the coherent channel noise) and step 3 (the correction angles).

A correction on the channel qubits q2, q3 is six numbers, the U(theta, phi, lambda) angles of each qubit, as in
qiskit's ``UGate``. Composition of a new rotation with an existing correction is done on SU(2) and converted back
with qiskit's ``OneQubitEulerDecomposer("U3")``, so the angles agree with what the hardware circuits were bound to.
"""
from __future__ import annotations

from functools import lru_cache
from typing import Mapping, Sequence

import numpy as np

from .paulis import I2, SIGMA, X, Y, Z

#: Stage 2 coherent channel noise: rotation axis and angle (rad) of each channel qubit
NOMINAL_NOISE: dict[int, dict] = {2: {"axis": (0.1, -0.5, -0.1), "angle": 0.5},
                                  3: {"axis": (-0.7, 1.0, 0.1), "angle": 1.3}}
#: 120-degree rotation about (1, 1, 1): x -> y -> z -> x
R_CYC: np.ndarray = np.array([[0.0, 0, 1], [1, 0, 0], [0, 1, 0]])


def rotation(axis: Sequence[float], angle: float) -> np.ndarray:
    """exp(-i angle/2 n.sigma): rotation of the Bloch sphere by ``angle`` about the unit axis n."""
    n = np.asarray(axis, float) / np.linalg.norm(axis)
    return np.cos(angle / 2) * I2 - 1j * np.sin(angle / 2) * sum(c * p for c, p in zip(n, SIGMA))


def axis_rotation(axis: int, angle: float) -> np.ndarray:
    """Rotation about x (0), y (1) or z (2)."""
    return np.cos(angle / 2) * I2 - 1j * np.sin(angle / 2) * SIGMA[axis]


def u_matrix(theta: float, phi: float, lam: float) -> np.ndarray:
    """Matrix of qiskit's U(theta, phi, lambda)."""
    c, s = np.cos(theta / 2), np.sin(theta / 2)
    return np.array([[c, -np.exp(1j * lam) * s], [np.exp(1j * phi) * s, np.exp(1j * (phi + lam)) * c]])


def so3(V: np.ndarray) -> np.ndarray:
    """Bloch rotation of a single-qubit unitary: V s_b V^dag = sum_a R[a, b] s_a."""
    Vd = V.conj().T
    return np.array([[0.5 * np.trace(SIGMA[a] @ V @ SIGMA[b] @ Vd).real for b in range(3)] for a in range(3)])


def su2(R: np.ndarray) -> np.ndarray:
    """SU(2) matrix w I - i (x X + y Y + z Z) of a rotation matrix (quaternion by Shepperd's method)."""
    t = np.trace(R)
    if t > 0:
        s = 2 * np.sqrt(t + 1)
        w, x, y, z = s / 4, (R[2, 1] - R[1, 2]) / s, (R[0, 2] - R[2, 0]) / s, (R[1, 0] - R[0, 1]) / s
    elif R[0, 0] > R[1, 1] and R[0, 0] > R[2, 2]:
        s = 2 * np.sqrt(1 + R[0, 0] - R[1, 1] - R[2, 2])
        w, x, y, z = (R[2, 1] - R[1, 2]) / s, s / 4, (R[0, 1] + R[1, 0]) / s, (R[0, 2] + R[2, 0]) / s
    elif R[1, 1] > R[2, 2]:
        s = 2 * np.sqrt(1 + R[1, 1] - R[0, 0] - R[2, 2])
        w, x, y, z = (R[0, 2] - R[2, 0]) / s, (R[0, 1] + R[1, 0]) / s, s / 4, (R[1, 2] + R[2, 1]) / s
    else:
        s = 2 * np.sqrt(1 + R[2, 2] - R[0, 0] - R[1, 1])
        w, x, y, z = (R[1, 0] - R[0, 1]) / s, (R[0, 2] + R[2, 0]) / s, (R[1, 2] + R[2, 1]) / s, s / 4
    return w * I2 - 1j * (x * X + y * Y + z * Z)


@lru_cache(maxsize=1)
def _u3_decomposer():
    from qiskit.synthesis import OneQubitEulerDecomposer
    return OneQubitEulerDecomposer("U3")


def angles_of(V: np.ndarray) -> np.ndarray:
    """U(theta, phi, lambda) angles of a single-qubit unitary (global phase dropped)."""
    return np.array(_u3_decomposer().angles(np.asarray(V, complex)))


def compose(params: Sequence[float] | None, R: Mapping[int, np.ndarray]) -> np.ndarray:
    """Angles of the new correction R[q] . V_old[q] on q2 and q3 (``params`` None = no correction yet)."""
    p = np.zeros(6) if params is None else np.asarray(params, float)
    return np.r_[angles_of(su2(R[2]) @ u_matrix(*p[:3])), angles_of(su2(R[3]) @ u_matrix(*p[3:6]))]


def correction_rotations(params: Sequence[float] | None) -> dict[int, np.ndarray]:
    """Bloch rotations of a six-angle correction (None = identity)."""
    if params is None:
        return {2: np.eye(3), 3: np.eye(3)}
    return {2: so3(u_matrix(*params[:3])), 3: so3(u_matrix(*params[3:6]))}


def noise_unitaries(kind: str = "N1") -> dict[int, np.ndarray]:
    """N0 = no injected noise; N1 = the Stage 2 coherent rotations (q2: 0.5 rad, q3: 1.3 rad)."""
    if kind == "N0":
        return {2: I2.copy(), 3: I2.copy()}
    if kind != "N1":
        raise KeyError(kind)
    return {q: rotation(p["axis"], p["angle"]) for q, p in NOMINAL_NOISE.items()}


def oracle_angles(kind: str = "N1") -> list[float]:
    """Correction angles of the exact inverse of the injected noise (what a perfect calibrator would find)."""
    U = noise_unitaries(kind)
    return [float(v) for v in angles_of(U[2].conj().T)] + [float(v) for v in angles_of(U[3].conj().T)]


def rotation_angle_axis(R: np.ndarray) -> tuple[float, list[float]]:
    """Angle (rad) and unit axis of a rotation matrix (axis (0, 0, 1) when the angle is 0)."""
    ang = float(np.arccos(np.clip((np.trace(R) - 1) / 2, -1, 1)))
    v = np.array([R[2, 1] - R[1, 2], R[0, 2] - R[2, 0], R[1, 0] - R[0, 1]])
    nv = np.linalg.norm(v)
    return ang, (v / nv).tolist() if nv > 1e-12 else [0, 0, 1]
