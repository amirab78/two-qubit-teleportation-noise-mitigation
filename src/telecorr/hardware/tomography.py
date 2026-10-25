"""Two-qubit state tomography of Bob's qubits from 9 Pauli-basis count tables.

Roadmap: step 8 (real hardware).

* ``linear``: every Pauli expectation <P4 P5> is averaged over the bases compatible with it, and
  rho = sum <P4 P5> sigma_P5 (x) sigma_P4 / 4. Its fidelity F = <psi|rho|psi> is *linear* in the 36 observed
  frequencies (unbiased); this is the main hardware estimate.
* ``mle``: the closest physical state, by projecting the eigenvalues of the linear estimate onto the probability
  simplex (Smolin, Gambetta and Smith, Phys. Rev. Lett. 108, 070502 (2012)). It is biased low for nearly pure states
  and is used for purity, Bloch vectors and correlators.
* error bars: parametric bootstrap, a multinomial resample of each basis at its observed frequencies.
"""
from __future__ import annotations

from typing import Mapping

import numpy as np

from ..paulis import KEYS4, PAULI, TOMO_BASES, correlation_matrix

Counts9 = Mapping[str, Mapping[str, float]]


def _sign(p4: str, p5: str, key: str) -> int:
    b5, b4 = int(key[0]), int(key[1])
    return (-1) ** ((b4 if p4 != "I" else 0) + (b5 if p5 != "I" else 0))


def linear(counts9: Counts9) -> tuple[np.ndarray, dict[str, float]]:
    """Linear-inversion state (qiskit order: acts on q5 (x) q4) and the 16 Pauli expectations."""
    exp = {}
    for p4 in "IXYZ":
        for p5 in "IXYZ":
            vals = []
            for bases, cnt in counts9.items():
                if (p4 != "I" and bases[0] != p4) or (p5 != "I" and bases[1] != p5):
                    continue
                n = sum(cnt.values())
                vals.append(sum(_sign(p4, p5, k) * v for k, v in cnt.items()) / n)
            exp[p4 + p5] = float(np.mean(vals))
    rho = sum(exp[p4 + p5] * np.kron(PAULI[p5], PAULI[p4]) for p4 in "IXYZ" for p5 in "IXYZ") / 4
    return rho, exp


def mle(counts9: Counts9) -> np.ndarray:
    """Physical state closest to the linear estimate (eigenvalue projection onto the simplex)."""
    rho, _ = linear(counts9)
    rho = (rho + rho.conj().T) / 2
    mu, V = np.linalg.eigh(rho)
    mu, V = mu[::-1], V[:, ::-1]
    lam = np.zeros(4); a = 0.0; i = 4
    while i > 0 and mu[i - 1] + a / i < 0:
        a += mu[i - 1]; i -= 1
    lam[:i] = mu[:i] + a / i if i > 0 else 0
    return (V * lam) @ V.conj().T


def fidelity(rho: np.ndarray, psi: np.ndarray) -> float:
    return float(np.real(np.conj(psi) @ rho @ psi))


def summary(rho: np.ndarray, psi: np.ndarray | None = None) -> dict:
    """Purity, (fidelity), Bloch vectors of q4 and q5 and the correlators <XX>, <YY>, <ZZ>."""
    out = {"purity": float(np.real(np.trace(rho @ rho)))}
    if psi is not None:
        out["F"] = fidelity(rho, psi)
    I2 = PAULI["I"]
    for name, embed in (("q4", lambda P: np.kron(I2, P)), ("q5", lambda P: np.kron(P, I2))):
        out[f"bloch_{name}"] = [float(np.real(np.trace(rho @ embed(PAULI[a])))) for a in "XYZ"]
    for a in "XYZ":
        out[f"<{a}{a}>"] = float(np.real(np.trace(rho @ np.kron(PAULI[a], PAULI[a]))))
    return out


def bootstrap_sigma(counts9: Counts9, psi: np.ndarray, n: int = 300, seed: int = 1) -> tuple[float, float]:
    """(sigma of the linear F, sigma of the MLE F) from ``n`` multinomial resamples of every basis.
    The random draws are made in a fixed order (resample, then basis in table order), so a seed fixes the result."""
    rng = np.random.default_rng(seed)
    arr = {bs: np.array([c.get(k, 0) for k in KEYS4], float) for bs, c in counts9.items()}
    Fl, Fm = [], []
    for _ in range(n):
        res = {}
        for bs, v in arr.items():
            p = np.clip(v, 0, None) / np.clip(v, 0, None).sum()
            res[bs] = dict(zip(KEYS4, rng.multinomial(int(round(v.sum())), p)))
        Fl.append(fidelity(linear(res)[0], psi)); Fm.append(fidelity(mle(res), psi))
    return float(np.std(Fl, ddof=1)), float(np.std(Fm, ddof=1))


def analyse(counts9: Counts9, psi: np.ndarray, boot: int = 300, seed: int = 1) -> dict:
    """Bob's state from 9 count tables: F (linear), F_mle, purity / Bloch / correlators of the MLE state, bootstrap
    sigmas (``boot`` resamples, 0 = none) and the MLE density matrix."""
    rho = mle(counts9)
    lin, _ = linear(counts9)
    out = summary(rho, psi)
    out["F_mle"] = out.pop("F")
    out["F"] = fidelity(lin, psi)
    if boot:
        out["F_sigma"], out["F_mle_sigma"] = bootstrap_sigma(counts9, psi, boot, seed)
    out["rho_real"] = np.real(rho).tolist(); out["rho_imag"] = np.imag(rho).tolist()
    return out


# ------------------------------------------------------------- the linear F as an explicit linear form
def linear_weights(psi: np.ndarray) -> dict[str, np.ndarray]:
    """w[B][k] such that F_lin = sum_B sum_k w[B][k] f[B][k] with f the observed frequencies of basis B."""
    c = correlation_matrix(psi)
    L = "IXYZ"
    w = {B: np.zeros(4) for B in TOMO_BASES}
    for i4, p4 in enumerate(L):
        for i5, p5 in enumerate(L):
            comp = [B for B in TOMO_BASES if (p4 == "I" or B[0] == p4) and (p5 == "I" or B[1] == p5)]
            for B in comp:
                for k, key in enumerate(KEYS4):
                    w[B][k] += c[i4, i5] * _sign(p4, p5, key) / len(comp) / 4
    return w


def linear_F(counts9: Counts9, w: Mapping[str, np.ndarray]) -> float:
    return float(sum(np.dot(w[B], [counts9[B].get(k, 0) for k in KEYS4]) / sum(counts9[B].values()) for B in TOMO_BASES))


def linear_F_variance(counts9: Counts9, w: Mapping[str, np.ndarray]) -> float:
    """Exact multinomial variance of F_lin at the plug-in frequencies (bases independent)."""
    v = 0.0
    for B in TOMO_BASES:
        n = np.array([counts9[B].get(k, 0) for k in KEYS4], float)
        N = n.sum(); f = n / N
        v += (np.dot(w[B] ** 2, f) - np.dot(w[B], f) ** 2) / N
    return float(v)
