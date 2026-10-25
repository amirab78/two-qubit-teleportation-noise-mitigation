"""Proven results behind the hardware analysis (Stage 9, Part 9.8), checked on the re-implementation."""
import itertools
import math

import numpy as np

from telecorr.hardware.analysis import fieller
from telecorr.hardware.frame import frame_bit, frame_keys
from telecorr.hardware.readout import confusion_from_counts
from telecorr.hardware.tomography import (linear, linear_F, linear_F_variance, linear_weights, mle)
from telecorr.paulis import KEYS4, PAULI, TOMO_BASES
from telecorr.states import H, S, haar_state, message_state

BASIS_U = {"Z": PAULI["I"], "X": H, "Y": H @ S.conj().T}            # rotation before a Z measurement


def projector(basis, b):
    return (PAULI["I"] + (-1) ** b * PAULI[basis]) / 2


def basis_probs(rho, bs):
    """Exact probabilities {"b5b4": p} of measuring Bob's (q4, q5) in the bases bs = (q4 basis, q5 basis)."""
    U = np.kron(BASIS_U[bs[1]], BASIS_U[bs[0]])
    p = np.real(np.diag(U @ rho @ U.conj().T))
    return {f"{k >> 1}{k & 1}": float(p[k]) for k in range(4)}


def test_T8_1_pauli_frame_equals_live_correction():
    """T8.1: for every basis B and correction Z^z X^x, (Z^z X^x)^dag Pi_b^B (Z^z X^x) = Pi_{b xor f}^B with
    f = frame_bit(0, B, z, x): measuring and flipping in software = correcting before the measurement."""
    for B, z, x, b in itertools.product("XYZ", (0, 1), (0, 1), (0, 1)):
        C = np.linalg.matrix_power(PAULI["Z"], z) @ np.linalg.matrix_power(PAULI["X"], x)
        lhs = C.conj().T @ projector(B, b) @ C
        assert np.allclose(lhs, projector(B, b ^ frame_bit(0, B, z, x)))


def test_T8_1_frame_keys_bit_layout():
    """Alice's bit 0 / 1 = Z corrections of q4 / q5, bits 2 / 3 = X corrections (Bob key = b5 * 2 + b4)."""
    for a in range(16):
        for bob in range(4):
            k = frame_keys(np.array([a]), np.array([bob]), "ZX")[0]
            # q4 measured in Z: flipped by its X bit (bit 2); q5 in X: flipped by its Z bit (bit 1)
            assert (k & 1) == ((bob & 1) ^ ((a >> 2) & 1))
            assert (k >> 1) == (((bob >> 1) & 1) ^ ((a >> 1) & 1))


def test_T8_2_readout_law():
    """T8.2: a bit with flip probabilities e0 = P(1|0), e1 = P(0|1) measures lambda_hat = (1 - e0 - e1) lambda + (e1 - e0)."""
    rng = np.random.default_rng(8002)
    for _ in range(200):
        p0, e0, e1 = rng.uniform(0, 1), rng.uniform(0, 0.2), rng.uniform(0, 0.2)
        lam = 2 * p0 - 1
        p0_hat = p0 * (1 - e0) + (1 - p0) * e1
        assert math.isclose(2 * p0_hat - 1, (1 - e0 - e1) * lam + (e1 - e0), abs_tol=1e-14)


def test_readout_calibration_counts():
    """confusion_from_counts recovers e0, e1 and Bob's matrix from pattern counts built with known errors."""
    e0 = [0.01, 0.02, 0.0, 0.03, 0.04, 0.05]; e1 = [0.02, 0.0, 0.01, 0.05, 0.03, 0.06]; N = 100000
    cc = {}
    for pat in ("000000", "111111", "010101", "101010"):
        reg = {}
        for name, qs in (("alice", range(4)), ("bob", (4, 5))):
            counts = {}
            for bits in itertools.product((0, 1), repeat=len(qs)):        # bits[j] = outcome of qubit qs[j]
                p = 1.0
                for j, q in enumerate(qs):
                    prep = int(pat[q]); flip = e0[q] if prep == 0 else e1[q]
                    p *= flip if bits[j] != prep else 1 - flip
                key = "".join(str(b) for b in reversed(bits))             # qiskit order: last qubit first
                counts[key] = p * N
            reg[name] = counts
        cc[pat] = reg
    r = confusion_from_counts(cc)
    assert np.allclose(r["e0"], e0) and np.allclose(r["e1"], e1)
    M = np.array(r["M_bob"])
    assert np.allclose(M.sum(0), 1)
    assert math.isclose(M[1, 0], e0[4] * (1 - e0[5]))                    # prepared 00, measured b4 = 1


def test_T8_3_linear_tomography_is_exact_on_probabilities():
    """T8.3: with exact probabilities the linear estimate returns the state, so F_lin is unbiased."""
    rng = np.random.default_rng(8003)
    for _ in range(20):
        psi = haar_state(rng)
        A = rng.normal(size=(4, 4)) + 1j * rng.normal(size=(4, 4))
        rho = A @ A.conj().T; rho /= np.trace(rho)
        c9 = {bs: basis_probs(rho, bs) for bs in TOMO_BASES}
        assert np.allclose(linear(c9)[0], rho, atol=1e-12)
        w = linear_weights(psi)
        assert math.isclose(linear_F(c9, w), float(np.real(psi.conj() @ rho @ psi)), abs_tol=1e-12)


def test_T8_3_linear_F_is_linear_in_frequencies_and_variance_exact():
    """F_lin = sum w f (the weights reproduce tomography's F on any counts); its multinomial variance formula agrees
    with a seeded Monte Carlo within 5 %."""
    rng = np.random.default_rng(8004)
    psi = message_state("bell_30_70")
    rho = 0.8 * np.outer(psi, psi.conj()) + 0.2 * np.eye(4) / 4
    w = linear_weights(psi)
    probs = {bs: basis_probs(rho, bs) for bs in TOMO_BASES}
    N = 3000
    sims = []
    for _ in range(3000):
        c9 = {bs: dict(zip(KEYS4, rng.multinomial(N, [max(probs[bs][k], 0) for k in KEYS4]))) for bs in TOMO_BASES}
        F_w = linear_F(c9, w)
        assert math.isclose(F_w, float(np.real(psi.conj() @ linear(c9)[0] @ psi)), abs_tol=1e-12)
        sims.append(F_w)
    c_exp = {bs: {k: probs[bs][k] * N for k in KEYS4} for bs in TOMO_BASES}
    v = linear_F_variance(c_exp, w)
    assert abs(np.var(sims, ddof=1) / v - 1) < 0.05
    assert abs(np.mean(sims) - float(np.real(psi.conj() @ rho @ psi))) < 4 * math.sqrt(v / len(sims))


def test_T8_3_mle_is_a_density_matrix_and_fixes_physical_states():
    rng = np.random.default_rng(8005)
    for _ in range(50):
        A = rng.normal(size=(4, 4)) + 1j * rng.normal(size=(4, 4))
        rho = A @ A.conj().T; rho /= np.trace(rho)
        c9 = {bs: basis_probs(rho, bs) for bs in TOMO_BASES}
        assert np.allclose(mle(c9), rho, atol=1e-12)                  # already physical: unchanged
        noisy = {bs: {k: max(v + rng.normal(0, 0.05), 0) for k, v in d.items()} for bs, d in c9.items()}
        r = mle(noisy)
        ev = np.linalg.eigvalsh(r)
        assert ev.min() > -1e-12 and math.isclose(np.trace(r).real, 1, abs_tol=1e-12)


def test_P8_10_delta_method_and_fieller():
    """P8.10: Var f = (s_c^2 + f^2 s_a^2 + (1 - f)^2 s_n^2) / D^2 for f = (c - n)/(a - n) (numeric derivative check),
    and the Fieller endpoints solve (N - f D)^2 = z^2 (Var N - 2 f Cov + f^2 Var D)."""
    rng = np.random.default_rng(8010)
    for _ in range(200):
        c, n, a = rng.uniform(0.75, 0.9), rng.uniform(0.4, 0.7), rng.uniform(0.8, 0.95)
        sc, sn, sa = rng.uniform(0.002, 0.01, 3)
        f = (c - n) / (a - n)
        h = 1e-7
        g = [((c + h - n) / (a - n) - f) / h, ((c - n - h) / (a - n - h) - f) / h, ((c - n) / (a + h - n) - f) / h]
        var_num = (g[0] * sc) ** 2 + (g[1] * sn) ** 2 + (g[2] * sa) ** 2
        var_formula = (sc ** 2 + f ** 2 * sa ** 2 + (1 - f) ** 2 * sn ** 2) / (a - n) ** 2
        assert math.isclose(var_num, var_formula, rel_tol=1e-4)
        Nv, Dv = c - n, a - n
        vN, vD, cND = sc ** 2 + sn ** 2, sa ** 2 + sn ** 2, sn ** 2
        lo, hi = fieller(Nv, Dv, vN, vD, cND, 1.96)
        assert lo < f < hi
        for r in (lo, hi):
            assert math.isclose((Nv - r * Dv) ** 2, 1.96 ** 2 * (vN - 2 * r * cND + r * r * vD), rel_tol=1e-9, abs_tol=1e-15)
