"""Correlated (crosstalk) noise between the two channel qubits: the two-qubit twirl, ceilings and TOMO-75.

Roadmap: step 6 (correlated noise).

Teleportation twirls any channel on (q2, q3) into a *two-qubit* Pauli channel. Write the unital 15 x 15 block of
the channel's Pauli transfer matrix as T (T[P, Q] = tr(s_P N(s_Q)) / 4 over the 15 non-identity Paulis). A
correction W (any two-qubit unitary, adjoint matrix Ad(W)) gives Bob the Pauli eigenvalues lambda_P = (Ad(W) T)_PP
(Theorem T6.2), so

    p_II = (1 + tr(Ad(W) T)) / 16,        F_avg = (4 p_II + 1) / 5                         (Theorem T6.3).

* A *local* correction (one unitary per qubit, what TOMO-9 learns) maximises tr(R2 T22) + tr(R3 T33) +
  tr((R2 (x) R3) Tcc): the correlated block couples the legs; alternating Procrustes solves it.
* The best *two-qubit* correction has no SVD closed form; p_II(W) <= lambda_max(chi) (chi = process matrix) is a
  certificate, the top eigenvector gives a start, and Newton steps on SU(4) polish it. When the incoherent part is
  symmetric positive semidefinite the inverse of the unitary part is optimal (Theorem T6.5).
* TOMO-75: a circuit with a Clifford setting S and a stabilizer probe reads exactly three entries of M = Ad(W) T;
  75 circuits read each of the 225 entries exactly once (product probes alone need >= 135; Theorem T6.11).

Pauli labels (j, k): s_j on q2 (message q0 -> Bob q4) and s_k on q3 (q1 -> q5); matrices in qiskit order q3 (x) q2.
"""
from __future__ import annotations

import json
from functools import lru_cache
from typing import Mapping, Sequence

import numpy as np
from scipy.linalg import expm

from .ceiling import procrustes
from .channels import incoherent_kraus, kraus_from_choi, choi, leg_parameters
from .paulis import I2, SIGMA
from .rotations import NOMINAL_NOISE, rotation, su2, so3, u_matrix
from .tomography import DESIGNS

P1 = (I2,) + SIGMA
LAB15 = [(j, 0) for j in (1, 2, 3)] + [(0, k) for k in (1, 2, 3)] + [(j, k) for j in (1, 2, 3) for k in (1, 2, 3)]
IDX = {lab: i for i, lab in enumerate(LAB15)}
S2, S3, SC = slice(0, 3), slice(3, 6), slice(6, 15)
NAMES = "IXYZ"
LAB_TXT = [NAMES[j] + NAMES[k] for j, k in LAB15]
LAB16 = [(j, k) for j in range(4) for k in range(4)]                   # natural order 4 j + k (0 = II)


def pauli2(j: int, k: int) -> np.ndarray:
    """s_j on q2 and s_k on q3, in qiskit order (q3 (x) q2)."""
    return np.kron(P1[k], P1[j])


SIG15 = np.stack([pauli2(j, k) for j, k in LAB15])
SIG16 = np.stack([pauli2(j, k) for j, k in LAB16])
_C1 = np.array([[1 if (a == 0 or b == 0 or a == b) else -1 for b in range(4)] for a in range(4)])
CHI16 = np.array([[_C1[j, l] * _C1[k, m] for (l, m) in LAB16] for (j, k) in LAB16])     # commutation signs


def ptm15(kraus: Sequence[np.ndarray]) -> np.ndarray:
    """Unital 15 x 15 block of the Pauli transfer matrix of a two-qubit channel."""
    K = np.asarray(kraus)
    NQ = np.einsum("iab,qbc,idc->qad", K, SIG15, K.conj())
    return np.einsum("pab,qba->pq", SIG15, NQ).real / 4


def adj(W: np.ndarray) -> np.ndarray:
    """Adjoint matrix of a two-qubit unitary in the 15-Pauli basis: W s_Q W^dag = sum_P Ad[P, Q] s_P."""
    return ptm15([W])


def _ad_generators() -> np.ndarray:
    out = np.zeros((15, 15, 15))
    for g in range(15):
        G = SIG15[g]
        comm = np.einsum("ab,qbc->qac", G, SIG15) - np.einsum("qab,bc->qac", SIG15, G)
        out[g] = np.einsum("pab,qba->pq", SIG15, -0.5j * comm).real / 4
    return out


AD = _ad_generators()
AD_SYM = 0.5 * (np.einsum("gpr,hrq->ghpq", AD, AD) + np.einsum("hpr,grq->ghpq", AD, AD))


def su4_step(y: np.ndarray) -> np.ndarray:
    """exp(-i/2 sum_G y_G s_G): local coordinates on SU(4)."""
    return expm(-0.5j * np.einsum("g,gab->ab", np.asarray(y, float), SIG15))


def haar_su4(rng: np.random.Generator) -> np.ndarray:
    Zm = (rng.normal(size=(4, 4)) + 1j * rng.normal(size=(4, 4))) / np.sqrt(2)
    Q, R = np.linalg.qr(Zm)
    Q = Q * (np.diag(R) / np.abs(np.diag(R)))
    return Q / np.linalg.det(Q) ** 0.25


# ---------------------------------------------------------------- the correlated channel
def crosstalk_unitary(zz: float, ex: float) -> np.ndarray:
    """exp(-i/2 [zz ZZ + ex (XX + YY)/2]): always-on ZZ coupling and exchange between the two channel qubits."""
    return expm(-0.5j * (zz * pauli2(3, 3) + ex * (pauli2(1, 1) + pauli2(2, 2)) / 2))


#: the nominal channel of step 6: Stage 2 rotations, crosstalk (ZZ 0.6, exchange 0.3), correlated dephasing 0.05, M 0.05
NOMINAL6 = dict(coherent=NOMINAL_NOISE, zz=0.6, ex=0.3, pc=0.05, family="M", p=0.05)


def unitary_part(spec: Mapping, xt: float = 1.0) -> np.ndarray:
    """Local rotations, then the crosstalk unitary scaled by xt."""
    c = spec.get("coherent")
    U = np.eye(4, dtype=complex) if c is None else np.kron(rotation(c[3]["axis"], c[3]["angle"]),
                                                          rotation(c[2]["axis"], c[2]["angle"]))
    return crosstalk_unitary(xt * spec.get("zz", 0.0), xt * spec.get("ex", 0.0)) @ U


def kraus6(spec: Mapping = NOMINAL6, xt: float = 1.0) -> list[np.ndarray]:
    """Kraus operators (4 x 4, order q3 (x) q2) of the whole channel: unitary part, then correlated dephasing
    {sqrt(1 - pc) II, sqrt(pc) ZZ}, then the local incoherent family on each qubit (p on q2, 1.5 p on q3)."""
    U = unitary_part(spec, xt)
    pc = spec.get("pc", 0.0)
    corr = [np.sqrt(1 - pc) * np.eye(4), np.sqrt(pc) * pauli2(3, 3)] if pc > 0 else [np.eye(4)]
    loc = {q: incoherent_kraus(**leg_parameters(spec.get("family"), spec.get("p", 0.0), q)) for q in (2, 3)}
    ks = [(np.kron(b, a) @ C @ U).astype(complex) for b in loc[3] for a in loc[2] for C in corr]
    return kraus_from_choi(choi(ks))


def T15(spec: Mapping = NOMINAL6, xt: float = 1.0) -> np.ndarray:
    return ptm15(kraus6(spec, xt))


# ---------------------------------------------------------------- exact model of Bob's channel
def lam16(W: np.ndarray, T: np.ndarray) -> np.ndarray:
    """Pauli eigenvalues of Bob's twirled channel (16, natural order, lambda_II = 1)."""
    d = np.diag(adj(W) @ T)
    out = np.ones(16)
    for i, (j, k) in enumerate(LAB15):
        out[4 * j + k] = d[i]
    return out


def p16(W: np.ndarray, T: np.ndarray) -> np.ndarray:
    """Pauli error probabilities p[4 j + k] of Bob's channel."""
    return CHI16 @ lam16(W, T) / 16


def favg(W: np.ndarray, T: np.ndarray) -> float:
    """Haar-average fidelity (4 p_II + 1) / 5, p_II = (1 + tr(Ad(W) T)) / 16."""
    return float(0.2 + 0.8 * (1 + np.trace(adj(W) @ T)) / 16)


def fidelity(psi: np.ndarray, W: np.ndarray, T: np.ndarray) -> float:
    """Exact fidelity of a pure message: sum_jk p_jk |<psi| s_j (x) s_k |psi>|^2."""
    w = np.array([abs(np.conj(psi) @ SIG16[i] @ psi) ** 2 for i in range(16)])
    return float(p16(W, T) @ w)


# ---------------------------------------------------------------- ceilings
def local_V(R2: np.ndarray, R3: np.ndarray) -> np.ndarray:
    return np.kron(su2(R3), su2(R2))


def per_leg_correction(T: np.ndarray) -> np.ndarray:
    """What TOMO-9 on each leg converges to: Procrustes of each leg's own 3 x 3 block, ignoring the coupling."""
    return local_V(procrustes(T[S2, S2]), procrustes(T[S3, S3]))


def best_local(T: np.ndarray, n_random: int = 6, seed: int = 0, iters: int = 200) -> np.ndarray:
    """Best local correction by alternating Procrustes (each step exact, the objective never decreases)."""
    Tc = T[SC, SC].reshape(3, 3, 3, 3)
    obj = lambda R2, R3: np.trace(R2 @ T[S2, S2]) + np.trace(R3 @ T[S3, S3]) + np.trace(np.kron(R2, R3) @ T[SC, SC])  # noqa: E731
    rng = np.random.default_rng(seed)
    starts = [(procrustes(T[S2, S2]), procrustes(T[S3, S3]))]
    starts += [tuple(so3(u_matrix(*rng.uniform(-np.pi, np.pi, 3))) for _ in range(2)) for _ in range(n_random)]
    best = None
    for R2, R3 in starts:
        f_old = -np.inf
        for _ in range(iters):
            R2 = procrustes(T[S2, S2] + np.einsum("bd,cdab->ca", R3, Tc))
            R3 = procrustes(T[S3, S3] + np.einsum("ac,cdab->db", R2, Tc))
            f = obj(R2, R3)
            if f - f_old < 1e-13:
                break
            f_old = f
        if best is None or f > best[0]:
            best = (f, R2, R3)
    return local_V(best[1], best[2])


@lru_cache(maxsize=1)
def _chi_inverse() -> np.ndarray:
    A = np.zeros((256, 256), complex)
    for P in range(16):
        for Q in range(16):
            A[16 * P + Q] = np.einsum("ab,mbc,cd,nda->mn", SIG16[P], SIG16, SIG16[Q],
                                      SIG16.conj().transpose(0, 2, 1)).ravel() / 4
    return np.linalg.inv(A)


def chi_of(T: np.ndarray) -> np.ndarray:
    """Process matrix chi (natural Pauli order) of the unital channel with Pauli-transfer block T."""
    R = np.zeros((16, 16)); R[0, 0] = 1
    perm = [4 * j + k for j, k in LAB15]
    R[np.ix_(perm, perm)] = T
    chi = (_chi_inverse() @ R.ravel()).reshape(16, 16)
    return (chi + chi.conj().T) / 2


def favg_bound(T: np.ndarray) -> float:
    """Upper bound on F_avg over ALL corrections: 0.2 + 0.8 lambda_max(chi)."""
    return float(0.2 + 0.8 * np.linalg.eigvalsh(chi_of(T))[-1])


def newton_su4(T: np.ndarray, V0: np.ndarray, iters: int = 60, tol: float = 1e-12) -> np.ndarray:
    """Maximise tr(Ad(V) T) over SU(4) by Newton steps in local coordinates (exact gradient and Hessian)."""
    V = np.asarray(V0, complex)
    f = lambda V: np.trace(adj(V) @ T)                       # noqa: E731
    fV = f(V)
    for _ in range(iters):
        Xm = adj(V) @ T
        g = np.einsum("gpq,qp->g", AD, Xm)
        if np.abs(g).max() < tol:
            break
        Hm = np.einsum("ghpq,qp->gh", AD_SYM, Xm)
        ev = np.linalg.eigvalsh(Hm)
        y = -np.linalg.solve(Hm, g) if ev[-1] < -1e-9 else 0.5 * g / max(np.abs(ev).max(), 1.0)
        step = 1.0
        while step > 1e-6:
            Vn = su4_step(step * y) @ V
            fn = f(Vn)
            if fn >= fV - 1e-14:
                V, fV = Vn, fn
                break
            step /= 2
        else:
            break
    return V


def best_su4(T: np.ndarray) -> np.ndarray:
    """Best two-qubit correction: start from the top eigenvector of chi (projected on the unitaries), then Newton."""
    w, v = np.linalg.eigh(chi_of(T))
    U, _, Wh = np.linalg.svd(np.einsum("m,mab->ab", v[:, -1].conj(), SIG16))
    return newton_su4(T, U @ Wh)


def ceilings(T: np.ndarray) -> dict:
    """Corrections of increasing power: none, per-leg Procrustes (TOMO-9), best local, best two-qubit."""
    return {"no correction": np.eye(4), "per-leg Procrustes (TOMO-9)": per_leg_correction(T),
            "best local": best_local(T), "best two-qubit": best_su4(T)}


def zz_local_ceiling(phi: float) -> float:
    """Pure ZZ crosstalk exp(-i phi/2 ZZ): the best LOCAL correction reaches p_II = max(cos^2, sin^2)(phi/2) (T6.7)."""
    return float(max(np.cos(phi / 2) ** 2, np.sin(phi / 2) ** 2))


# ---------------------------------------------------------------- two-qubit Cliffords, probe groups and TOMO-75
def _gate_mats() -> dict[str, np.ndarray]:
    Hm = np.array([[1, 1], [1, -1]]) / np.sqrt(2); Sm = np.diag([1, 1j])
    cx = np.zeros((4, 4))                                    # control q2 (bit 0), target q3 (bit 1)
    for b3 in (0, 1):
        for b2 in (0, 1):
            cx[2 * (b3 ^ b2) + b2, 2 * b3 + b2] = 1
    return {"h0": np.kron(I2, Hm), "h1": np.kron(Hm, I2), "s0": np.kron(I2, Sm), "s1": np.kron(Sm, I2), "cx": cx}


def _perm_key(A: np.ndarray) -> tuple:
    return tuple(int(np.argmax(np.abs(A[:, q]))) for q in range(15))


@lru_cache(maxsize=1)
def clifford_table() -> tuple:
    """One two-qubit Clifford per symplectic class (720), found by breadth-first search over H, S and CX,
    sorted by its Pauli permutation."""
    G = _gate_mats()
    start = (np.eye(4, dtype=complex), [])
    seen = {_perm_key(adj(start[0])): start}
    frontier = [start]
    while frontier:
        nxt = []
        for U, gl in frontier:
            for name, M in G.items():
                Un = M @ U
                k = _perm_key(adj(Un))
                if k not in seen:
                    seen[k] = (Un, gl + [name]); nxt.append(seen[k])
        frontier = nxt
    return tuple(dict(U=seen[k][0], gates=seen[k][1], A=np.round(adj(seen[k][0])).astype(int)) for k in sorted(seen))


def _lab_product(p, q):
    sym = {0: (0, 0), 1: (1, 0), 2: (1, 1), 3: (0, 1)}
    back = {v: k for k, v in sym.items()}
    f = lambda a, b: back[(sym[a][0] ^ sym[b][0], sym[a][1] ^ sym[b][1])]     # noqa: E731
    return (f(p[0], q[0]), f(p[1], q[1]))


@lru_cache(maxsize=1)
def probe_groups() -> tuple:
    """The 15 two-qubit stabilizer groups {P1, P2, P1 P2}: 9 product and 6 entangled, each with a probe Clifford."""
    zq2, zq3 = IDX[(3, 0)], IDX[(0, 3)]
    out, seen = [], set()
    for C in clifford_table():
        g1 = int(np.argmax(np.abs(C["A"][:, zq2]))); g2 = int(np.argmax(np.abs(C["A"][:, zq3])))
        g3 = IDX[_lab_product(LAB15[g1], LAB15[g2])]
        key = frozenset((g1, g2, g3))
        if key not in seen:
            seen.add(key)
            prod = any(LAB15[g][1] == 0 for g in key) and any(LAB15[g][0] == 0 for g in key)
            out.append(dict(members=(g1, g2, g3), product=prod, C=C["U"],
                            text="{" + ", ".join(LAB_TXT[g] for g in (g1, g2, g3)) + "}"))
    out.sort(key=lambda d: (not d["product"], sorted(d["members"])))
    return tuple(out)


def circuit_rows(setting: int, group: int) -> list[tuple[int, int, int]]:
    """The three entries (P, Q, sign) of M = Ad(W) T that circuit (Clifford setting, probe group) reads."""
    A = clifford_table()[setting]["A"]
    rows = []
    for P in probe_groups()[group]["members"]:
        Q = int(np.argmax(np.abs(A[P]))); rows.append((P, Q, int(A[P, Q])))
    return rows


def design_matrix(design: Sequence[tuple[int, int]]) -> np.ndarray:
    """Linear map vec(M) -> measured lambdas (3 per circuit); vec index 15 Q + P <-> M[Q, P]."""
    A = np.zeros((3 * len(design), 225))
    for c, (s, g) in enumerate(design):
        for r, (P, Q, sg) in enumerate(circuit_rows(s, g)):
            A[3 * c + r, 15 * Q + P] = sg
    return A


def load_tomo75(path: str) -> list[tuple[int, int]]:
    """The TOMO-75 design saved by the study (75 pairs (Clifford setting, probe group))."""
    with open(path, encoding="utf8") as fh:
        return [tuple(x) for x in json.load(fh)["design"]]


def circuit_lambdas(design, M: np.ndarray) -> np.ndarray:
    return np.array([[sg * M[Q, P] for P, Q, sg in circuit_rows(s, g)] for s, g in design])


def outcome_probs(lam: Sequence[float]) -> np.ndarray:
    """P of Bob's keys "00", "01", "10", "11" (b5 b4) from the circuit's (lambda_g1, lambda_g2, lambda_g3)."""
    l1, l2, l3 = lam
    p = lambda b0, b1: (1 + (-1) ** b0 * l1 + (-1) ** b1 * l2 + (-1) ** (b0 + b1) * l3) / 4      # noqa: E731
    return np.clip(np.array([p(0, 0), p(1, 0), p(0, 1), p(1, 1)]), 0, None)


def sample_lambdas(rng: np.random.Generator, design, M: np.ndarray, N: int) -> np.ndarray:
    """Shot-noise estimates of the circuit lambdas with N shots per circuit (multinomial counts)."""
    out = []
    for lam in circuit_lambdas(design, M):
        p = outcome_probs(lam)
        n00, n01, n10, n11 = rng.multinomial(N, p / p.sum())
        out.append([(n00 + n10 - n01 - n11) / N, (n00 + n01 - n10 - n11) / N, (n00 + n11 - n01 - n10) / N])
    return np.array(out)


def estimate_M(design, lam_hat: np.ndarray) -> np.ndarray:
    """Least-squares estimate of M (15 x 15) from the measured lambdas (an exact inversion for TOMO-75)."""
    m = np.linalg.lstsq(design_matrix(design), np.asarray(lam_hat).ravel(), rcond=None)[0]
    return m.reshape(15, 15)


# ---------------------------------------------------------------- the crosstalk detector on TOMO-9 circuits
def tomo9_probs(T: np.ndarray, W: np.ndarray | None = None) -> list[np.ndarray]:
    """Outcome probabilities of the 9 TOMO-9 circuits (settings I, C, C^2 on both legs; probes Z, X, Y) on a
    correlated channel, after an optional current correction W."""
    W = np.eye(4) if W is None else W
    out = []
    for R in DESIGNS["tomo9"]:
        M = adj(local_V(R, R) @ W) @ T
        for a in (2, 0, 1):                                       # probe order |00> (z), |++> (x), |+i+i> (y)
            l1, l2, l3 = (M[IDX[(a + 1, 0)], IDX[(a + 1, 0)]], M[IDX[(0, a + 1)], IDX[(0, a + 1)]],
                          M[IDX[(a + 1, a + 1)], IDX[(a + 1, a + 1)]])
            out.append(outcome_probs((l1, l2, l3)))
    return out


def detector_signal(T: np.ndarray) -> np.ndarray:
    """D = lambda_aa - lambda_a0 lambda_0a in the 9 TOMO-9 circuits: exactly 0 for any product channel."""
    out = []
    for R in DESIGNS["tomo9"]:
        M = adj(local_V(R, R)) @ T
        for a in range(3):
            l1, l2, l3 = M[IDX[(a + 1, 0)], IDX[(a + 1, 0)]], M[IDX[(0, a + 1)], IDX[(0, a + 1)]], M[IDX[(a + 1, a + 1)], IDX[(a + 1, a + 1)]]
            out.append(l3 - l1 * l2)
    return np.array(out)
