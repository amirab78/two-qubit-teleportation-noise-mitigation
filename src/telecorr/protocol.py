"""The six-qubit two-qubit teleportation protocol: exact simulation and the qiskit circuits.

Roadmap: step 1 (teleportation); every later step runs through the same protocol.

Qubits: q0, q1 = Alice's message; (q2, q4) and (q3, q5) = the two Bell pairs; q4, q5 = Bob. The channel acts on
Alice's halves q2, q3 after the pairs are shared, a correction box (one unitary per channel qubit) follows it, then
Alice measures (q0, q2) and (q1, q3) in the Bell basis and Bob applies

    q4:  X^{m2} then Z^{m0},        q5:  X^{m3} then Z^{m1}            (m_k = Alice's bit of qubit k).

:func:`teleport` is an exact density-matrix simulation of this protocol (64 x 64), branch by branch over Alice's 16
outcomes. It needs no simulator and is the reference for every model in the package. :func:`circuit` builds the
same protocol as a qiskit circuit (deferred-measurement or mid-circuit measurement with feed-forward) for Aer.
"""
from __future__ import annotations

from typing import Callable, Mapping, Sequence

import numpy as np

from .paulis import I2, X, Z
from .states import H, density

N_QUBITS = 6
DIM = 2 ** N_QUBITS
NOISE_LABEL = "channel noise"

#: Bob's Pauli frame, i.e. which of Alice's bits drive which correction. Wrong wirings are used in notebook 1.
FRAMES: dict[str, Callable[[Sequence[int]], tuple[np.ndarray, np.ndarray]]] = {
    "correct": lambda m: (_pow(Z, m[0]) @ _pow(X, m[2]), _pow(Z, m[1]) @ _pow(X, m[3])),
    "no correction": lambda m: (I2, I2),
    "X and Z bits exchanged": lambda m: (_pow(Z, m[2]) @ _pow(X, m[0]), _pow(Z, m[3]) @ _pow(X, m[1])),
    "only X corrections": lambda m: (_pow(X, m[2]), _pow(X, m[3])),
    "only Z corrections": lambda m: (_pow(Z, m[0]), _pow(Z, m[1])),
}


def _pow(P: np.ndarray, k: int) -> np.ndarray:
    return P if k else I2


def embed(U: np.ndarray, low: int, n: int = N_QUBITS) -> np.ndarray:
    """An operator on the consecutive qubits low, low + 1, ... (qiskit order inside U) as a 2^n matrix."""
    k = int(round(np.log2(U.shape[0])))
    return np.kron(np.kron(np.eye(2 ** (n - low - k)), U), np.eye(2 ** low))


def cx(control: int, target: int, n: int = N_QUBITS) -> np.ndarray:
    """CNOT as a 2^n permutation matrix."""
    M = np.zeros((2 ** n, 2 ** n))
    for i in range(2 ** n):
        M[i ^ (1 << target) if (i >> control) & 1 else i, i] = 1
    return M


_BELL_PAIRS = cx(3, 5) @ embed(H, 3) @ cx(2, 4) @ embed(H, 2)               # H q2, CX q2->q4, H q3, CX q3->q5
_ALICE_BASIS = embed(H, 1) @ cx(1, 3) @ embed(H, 0) @ cx(0, 2)              # CX q0->q2, H q0, CX q1->q3, H q1


def _channel(rho: np.ndarray, kraus: Sequence[np.ndarray], low: int) -> np.ndarray:
    Ks = [embed(np.asarray(K, complex), low) for K in kraus]
    return sum(K @ rho @ K.conj().T for K in Ks)


def _branches(rho_m: np.ndarray, legs, pair, correction, bob_box, frame: str) -> list[np.ndarray]:
    """Bob's corrected, UNnormalised state for each of Alice's 16 outcomes (linear in the message rho_m)."""
    zero = np.zeros((16, 16)); zero[0, 0] = 1
    rho = np.kron(zero, rho_m)                                    # q5..q2 in |0>, message on q1 q0
    rho = _BELL_PAIRS @ rho @ _BELL_PAIRS.T
    if legs is not None:
        rho = _channel(_channel(rho, legs[2], 2), legs[3], 3)
    if pair is not None:
        rho = _channel(rho, pair, 2)
    if correction is not None:
        W = correction if isinstance(correction, np.ndarray) else np.kron(correction[3], correction[2])
        rho = _channel(rho, [W], 2)
    rho = _ALICE_BASIS @ rho @ _ALICE_BASIS.T
    r4 = rho.reshape(4, 16, 4, 16)                                 # (Bob q5 q4, Alice q3..q0) for rows and columns
    B = np.eye(4) if bob_box is None else np.kron(bob_box[5], bob_box[4])
    out = []
    for m in range(16):
        c4, c5 = FRAMES[frame]([(m >> k) & 1 for k in range(4)])
        C = B @ np.kron(c5, c4)
        out.append(C @ r4[:, m, :, m] @ C.conj().T)
    return out


def teleport(message: np.ndarray, legs: Mapping[int, Sequence[np.ndarray]] | None = None,
             pair: Sequence[np.ndarray] | None = None, correction: Mapping[int, np.ndarray] | np.ndarray | None = None,
             bob_box: Mapping[int, np.ndarray] | None = None, frame: str = "correct") -> dict:
    """Exact teleportation of a two-qubit message (state vector or 4 x 4 density matrix).

    legs        Kraus operators of the channel on each channel qubit, ``{2: [...], 3: [...]}``;
    pair        or Kraus operators (4 x 4, order q3 (x) q2) of a two-qubit channel on (q2, q3);
    correction  the correction box after the channel: ``{2: V2, 3: V3}`` or one 4 x 4 unitary on (q2, q3);
    bob_box     a fixed unitary per Bob qubit applied after Bob's corrections, ``{4: B4, 5: B5}`` (notebook 2);
    frame       which Pauli frame Bob uses (:data:`FRAMES`).

    Returns ``{"p": Alice's 16 outcome probabilities (index m0 + 2 m1 + 4 m2 + 8 m3), "branches": Bob's normalised
    state for each outcome after his corrections, "bob": Bob's average state}``.
    """
    msg = np.asarray(message, complex)
    sub = _branches(density(msg) if msg.ndim == 1 else msg, legs, pair, correction, bob_box, frame)
    p = np.array([float(np.real(np.trace(b))) for b in sub])
    branches = [b / pm if pm > 1e-15 else np.zeros((4, 4), complex) for b, pm in zip(sub, p)]
    return {"p": p, "branches": branches, "bob": sum(sub)}


def choi(**kw) -> np.ndarray:
    """Choi matrix (16 x 16, trace 1) of the whole protocol as a channel message (q1 q0) -> Bob (q5 q4).
    Takes the keyword arguments of :func:`teleport` (without the message)."""
    J = np.zeros((16, 16), complex)
    for i in range(4):
        for j in range(4):
            E = np.zeros((4, 4), complex); E[i, j] = 1
            J += np.kron(E, sum(_branches(E, kw.get("legs"), kw.get("pair"), kw.get("correction"),
                                          kw.get("bob_box"), kw.get("frame", "correct")))) / 4
    return J


def average_fidelity(**kw) -> float:
    """Exact Haar-average fidelity of the protocol over all two-qubit messages: (4 F_e + 1) / 5 with the
    entanglement fidelity F_e = <Omega| J |Omega> of the Choi matrix (keyword arguments of :func:`teleport`)."""
    omega = np.eye(4).reshape(16) / 2
    Fe = float(np.real(omega @ choi(**kw) @ omega))
    return (4 * Fe + 1) / 5


def fidelity(target: np.ndarray, rho: np.ndarray) -> float:
    """Fidelity of a state rho with a target (state vector: <psi|rho|psi>; density matrix: Uhlmann)."""
    t = np.asarray(target, complex)
    if t.ndim == 1:
        return float(np.real(t.conj() @ rho @ t))
    from scipy.linalg import sqrtm
    s = sqrtm(t)
    return float(np.real(np.trace(sqrtm(s @ rho @ s))) ** 2)


def purity(rho: np.ndarray) -> float:
    return float(np.real(np.trace(rho @ rho)))


# ---------------------------------------------------------------- the same protocol as qiskit circuits
def prep_unitary(psi: np.ndarray) -> np.ndarray:
    """A unitary whose first column is psi (prepares psi from |00>)."""
    psi = np.asarray(psi, complex)
    M = np.eye(4, dtype=complex); M[:, 0] = psi
    Q, R = np.linalg.qr(M)
    return Q * (R[0, 0] / abs(R[0, 0])).conj()


def _prepare(qc, name_or_psi, a: int, b: int, inverse: bool = False) -> None:
    from qiskit.circuit.library import UnitaryGate
    from .states import THETA_BELL
    if isinstance(name_or_psi, str) and name_or_psi == "product_HT":
        for q in (a, b):
            if inverse:
                qc.tdg(q); qc.h(q)
            else:
                qc.h(q); qc.t(q)
    elif isinstance(name_or_psi, str) and name_or_psi == "bell_30_70":
        if inverse:
            qc.cx(a, b); qc.ry(-THETA_BELL, a)
        else:
            qc.ry(THETA_BELL, a); qc.cx(a, b)
    elif isinstance(name_or_psi, str) and name_or_psi in ("|00>", "|++>", "|+i+i>"):
        # probes are prepared qubit by qubit, so each of Bob's bits reads its own leg after un-preparation
        for q in (a, b):
            if name_or_psi == "|++>":
                qc.h(q)
            elif name_or_psi == "|+i+i>":
                if inverse:
                    qc.sdg(q); qc.h(q)
                else:
                    qc.h(q); qc.s(q)
    else:
        from .states import message_state
        psi = message_state(name_or_psi) if isinstance(name_or_psi, str) else name_or_psi
        U = prep_unitary(psi)
        qc.append(UnitaryGate(U.conj().T if inverse else U, label="unprep" if inverse else "prep"), [a, b])


def circuit(message, version: str = "measured", basis: str = "ZZ", unprepare: bool = False,
            channel: bool = True, correction: Sequence[float] | None = None):
    """The protocol as a qiskit circuit.

    version    ``"unitary"``: deferred measurement (Bob's corrections as CX / CZ), no measurement of Alice;
               ``"measured"``: Alice measured mid-circuit, Bob's corrections by classical feed-forward.
    basis      Bob's measurement basis (q4 basis, q5 basis), e.g. ``"ZZ"`` or ``"XY"``; ignored with unprepare.
    unprepare  apply the inverse of the message preparation to Bob's qubits and measure: P(00) = F (pure message).
    channel    insert the black-box identity gate labelled "channel noise" on (q2, q3) (noise added by the model).
    correction six U(theta, phi, lambda) angles on q2 and q3 after the channel (None = no correction box).

    The unitary version without measurement is meant for density-matrix simulation (``save_density_matrix``).
    """
    from qiskit import ClassicalRegister, QuantumCircuit, QuantumRegister
    from qiskit.circuit.library import UnitaryGate
    q = QuantumRegister(6, "q")
    alice, bob = ClassicalRegister(4, "alice"), ClassicalRegister(2, "bob")
    measured = version == "measured"
    qc = QuantumCircuit(q, alice, bob) if measured else QuantumCircuit(q, bob)
    _prepare(qc, message, 0, 1); qc.barrier()
    qc.h(2); qc.cx(2, 4); qc.h(3); qc.cx(3, 5); qc.barrier()
    if channel:
        qc.append(UnitaryGate(np.eye(4), label=NOISE_LABEL), [2, 3]); qc.barrier()
    if correction is not None:
        qc.u(*correction[:3], 2); qc.u(*correction[3:6], 3); qc.barrier()
    qc.cx(0, 2); qc.h(0); qc.cx(1, 3); qc.h(1); qc.barrier()
    if measured:
        qc.measure([0, 1, 2, 3], alice); qc.barrier()
        with qc.if_test((alice[2], 1)):
            qc.x(4)
        with qc.if_test((alice[0], 1)):
            qc.z(4)
        with qc.if_test((alice[3], 1)):
            qc.x(5)
        with qc.if_test((alice[1], 1)):
            qc.z(5)
    else:
        qc.cx(2, 4); qc.cz(0, 4); qc.cx(3, 5); qc.cz(1, 5)
    qc.barrier()
    if unprepare:
        _prepare(qc, message, 4, 5, inverse=True)
    else:
        for qb, b in zip((4, 5), basis):
            if b == "X":
                qc.h(qb)
            elif b == "Y":
                qc.sdg(qb); qc.h(qb)
    return qc


def noise_model(pair_kraus: Sequence[np.ndarray] | None):
    """Aer noise model that turns the "channel noise" box into the given two-qubit channel (order q3 (x) q2)."""
    from qiskit_aer.noise import NoiseModel, kraus_error
    nm = NoiseModel(basis_gates=["unitary", "u", "h", "x", "z", "cx", "cz", "ry", "t", "tdg", "s", "sdg",
                                 "measure", "barrier"])
    if pair_kraus is not None:
        nm.add_quantum_error(kraus_error([np.asarray(K, complex) for K in pair_kraus]), NOISE_LABEL, [2, 3])
    return nm


def bob_counts(qc, shots: int, seed: int, pair_kraus: Sequence[np.ndarray] | None = None) -> dict[str, int]:
    """Measure Bob's register with Aer (``"b5b4"`` keys); with ``pair_kraus`` the channel box is noisy.
    Different ``seed`` values give independent samples."""
    from qiskit import ClassicalRegister
    from qiskit_aer import AerSimulator
    c = qc.copy()
    if not c.clbits or not any(r.name == "bob" for r in c.cregs):
        c.add_register(ClassicalRegister(2, "bob"))
    bob = next(r for r in c.cregs if r.name == "bob")
    c.measure([4, 5], bob)
    # Aer seeds shot k with seed + k, so nearby seeds would share most of their shots: spread the seed first
    aer_seed = int(np.random.default_rng(seed).integers(2 ** 31))
    raw = AerSimulator(noise_model=noise_model(pair_kraus)).run(c, shots=shots, seed_simulator=aer_seed).result().get_counts()
    out: dict[str, int] = {}
    for key, n in raw.items():
        k = key.split()[0]                       # the bob register is printed first
        out[k] = out.get(k, 0) + n
    return out


def aer_bob_state(qc, pair_kraus: Sequence[np.ndarray] | None = None) -> np.ndarray:
    """Bob's density matrix from Aer's density-matrix simulator (unitary version of :func:`circuit`)."""
    from qiskit.quantum_info import partial_trace
    from qiskit_aer import AerSimulator
    c = qc.copy(); c.save_density_matrix()
    rho = AerSimulator(method="density_matrix", noise_model=noise_model(pair_kraus)).run(c).result().data(0)["density_matrix"]
    return np.asarray(partial_trace(rho, [0, 1, 2, 3]).data)
