"""The two test messages, the three calibration probes and random states.

Roadmap: step 1 (teleportation); the messages are used in every step.

* ``product_HT``: H then T on each message qubit, (|0> + e^{i pi/4}|1>)/sqrt(2) on q0 and on q1;
* ``bell_30_70``: Ry(theta) on q0 then CX q0 -> q1, sqrt(0.3)|00> + sqrt(0.7)|11>, theta = 2 arcsin(sqrt(0.7)).

State vectors are in qiskit order |q1 q0> (Bob receives them on |q5 q4>).
"""
from __future__ import annotations

import numpy as np

from .paulis import I2

H: np.ndarray = np.array([[1, 1], [1, -1]], dtype=complex) / np.sqrt(2)
S: np.ndarray = np.diag([1, 1j]).astype(complex)
T_GATE: np.ndarray = np.diag([1, np.exp(1j * np.pi / 4)]).astype(complex)
THETA_BELL: float = float(2 * np.arcsin(np.sqrt(0.7)))

MESSAGES: tuple[str, str] = ("product_HT", "bell_30_70")
#: one-qubit preparations from |0> used by the hardware circuits (H first)
PREP_U: dict[str, np.ndarray] = {"product_HT": T_GATE @ H, "0": I2, "+": H, "+i": S @ H}
#: the three counts-only calibration probes and the Bloch axis each one tests (same state on both qubits)
PROBES: tuple[tuple[str, str], ...] = (("|00>", "Z"), ("|++>", "X"), ("|+i+i>", "Y"))
PROBE_PREP: dict[str, str] = {"|00>": "0", "|++>": "+", "|+i+i>": "+i"}


def ry(theta: float) -> np.ndarray:
    c, s = np.cos(theta / 2), np.sin(theta / 2)
    return np.array([[c, -s], [s, c]], dtype=complex)


def message_state(name: str) -> np.ndarray:
    """State vector (qiskit order |q1 q0>) of a test message or a probe."""
    zero = np.array([1, 0], dtype=complex)
    if name == "product_HT":
        a = T_GATE @ H @ zero
        return np.kron(a, a)
    if name == "bell_30_70":
        a = ry(THETA_BELL) @ zero                    # q0
        psi = np.kron(np.array([1, 0], dtype=complex), a)   # q1 = |0>
        cx = np.eye(4, dtype=complex)[[0, 3, 2, 1]]  # CX q0 -> q1 in qiskit order: |01> <-> |11>
        return cx @ psi
    if name in PROBE_PREP:
        a = PREP_U[PROBE_PREP[name]] @ zero
        return np.kron(a, a)
    raise KeyError(name)


def density(psi: np.ndarray) -> np.ndarray:
    return np.outer(psi, np.conj(psi))


def haar_state(rng: np.random.Generator, dim: int = 4) -> np.ndarray:
    """A Haar-random pure state of dimension ``dim``."""
    v = rng.normal(size=dim) + 1j * rng.normal(size=dim)
    return v / np.linalg.norm(v)


def random_message(seed: int) -> np.ndarray:
    """The Haar-random message U|00> of qiskit's ``random_unitary(4, seed)`` (the random test messages of the study)."""
    from qiskit.quantum_info import random_unitary
    return np.asarray(random_unitary(4, seed=seed).data)[:, 0]
