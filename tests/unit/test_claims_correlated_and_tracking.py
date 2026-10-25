"""Proven results of steps 6-7 (Parts 9.6-9.7): two-qubit twirl, ceilings, TOMO-75, TOMO-9 and the Kalman tracker."""
import numpy as np
import pytest

from telecorr import correlated as X
from telecorr import protocol as P
from telecorr import tracking as Tr
from telecorr.ceiling import procrustes_batch
from telecorr.paulis import Z
from telecorr.reference import path
from telecorr.states import haar_state
from telecorr.tomography import DESIGNS, ls_T

RNG = np.random.default_rng(96)


def test_T6_2_T6_3_two_qubit_twirl_and_fidelities():
    for spec in (X.NOMINAL6, dict(X.NOMINAL6, zz=-0.8, ex=0.5, family="A", p=0.2)):
        Kr = X.kraus6(spec); T = X.ptm15(Kr)
        W = X.haar_su4(RNG)
        psi = haar_state(RNG)
        exact = P.fidelity(psi, P.teleport(psi, pair=Kr, correction=W)["bob"])
        assert exact == pytest.approx(X.fidelity(psi, W, T), abs=1e-12)
        assert P.average_fidelity(pair=Kr, correction=W) == pytest.approx(X.favg(W, T), abs=1e-12)


def test_T6_5_inverse_is_optimal_for_symmetric_psd_noise():
    for _ in range(5):
        lam = RNG.uniform(0.05, 1, 15)                         # a two-qubit Pauli channel with all lambda > 0
        U = X.haar_su4(RNG)
        T = np.diag(lam) @ X.adj(U)
        Uinv = U.conj().T
        assert X.favg(X.best_su4(T), T) <= X.favg(Uinv, T) + 1e-9
        assert all(X.favg(X.haar_su4(RNG), T) < X.favg(Uinv, T) for _ in range(50))


def test_T6_7_pure_zz_crosstalk_local_ceiling():
    for phi in np.linspace(0.1, np.pi, 7):
        T = X.ptm15([X.crosstalk_unitary(phi, 0.0)])
        pII = (1 + np.trace(X.adj(X.best_local(T)) @ T)) / 16
        assert pII == pytest.approx(X.zz_local_ceiling(phi), abs=1e-9)
        assert X.favg(X.crosstalk_unitary(phi, 0.0).conj().T, T) == pytest.approx(1, abs=1e-12)


def test_T6_11_tomo75_reads_every_entry_exactly_once():
    design = X.load_tomo75(path("6a_tomo75_design.json"))
    groups = X.probe_groups()
    assert len(X.clifford_table()) == 720 and len(groups) == 15 and sum(g["product"] for g in groups) == 9
    A = X.design_matrix(design)
    assert len(design) == 75 and sum(groups[g]["product"] for _, g in design) == 45
    assert np.allclose(A.T @ A, np.eye(225))                 # a signed permutation: exact and minimal
    M = X.T15() @ X.adj(X.haar_su4(RNG))
    assert np.abs(X.estimate_M(design, X.circuit_lambdas(design, M)) - M).max() < 1e-12


def test_P6_13_detector_signal_is_zero_for_product_channels():
    spec = dict(X.NOMINAL6, zz=0.0, ex=0.0, pc=0.0, family="A", p=0.2)
    assert np.abs(X.detector_signal(X.T15(spec))).max() < 1e-12
    assert np.abs(X.detector_signal(X.T15())).max() > 0.1


def test_T7_2_tomo9_reads_every_entry_of_T_exactly_once():
    C = DESIGNS["tomo9"]
    assert np.allclose(np.linalg.matrix_power(C[1], 3), np.eye(3))
    assert np.allclose(Tr.A9.T @ Tr.A9, np.eye(9)) and set(np.unique(Tr.A9)) == {0.0, 1.0}
    T = RNG.normal(size=(4, 3, 3))
    lam = np.einsum("jak,nka->nja", C, T)
    assert np.abs(ls_T(np.broadcast_to(C, (4, 3, 3, 3)), lam, ridge=0.0) - T).max() < 1e-12


def test_T7_3_procrustes_batch_equals_the_closed_form_value():
    T = RNG.normal(size=(200, 3, 3))
    R = procrustes_batch(T)
    s = np.linalg.svd(T, compute_uv=False)
    val = s[:, 0] + s[:, 1] + np.sign(np.linalg.det(T)) * s[:, 2]
    assert np.allclose(Tr.trace_RT(R, T), val, atol=1e-12)


def test_P7_8_sequential_kalman_updates_equal_the_batch_information_form():
    kf = Tr.KalmanTracker(4)
    T = RNG.normal(0, 0.5, (4, 3, 3)); N = np.full(4, 1000)
    kf.init_tomo9(np.stack([Tr.measure_lambda(RNG, np.broadcast_to(S, (4, 3, 3)), T, N) for S in DESIGNS["tomo9"]], 1).reshape(4, 9), N)
    S = np.broadcast_to(DESIGNS["tomo12"][1], (4, 3, 3))
    x0, P0 = kf.x.copy(), kf.P.copy()
    lam = Tr.measure_lambda(RNG, S, T, N)
    G = Tr.round_G(S)
    kf.update(S, lam, N)
    for m in range(4):                                         # batch form with the same plugged-in variances
        x, Pm = x0[m].copy(), P0[m].copy()
        for a in range(3):
            g = G[m, a]; pred = g @ x
            v = max(1 - np.clip(pred, -1, 1) ** 2, 1 / N[m]) / N[m]
            k = Pm @ g / (g @ Pm @ g + v)
            x = x + k * (lam[m, a] - pred); Pm = Pm - np.outer(k, g @ Pm)
        assert np.allclose(kf.x[m], x, atol=1e-10) and np.allclose(kf.P[m], Pm, atol=1e-10)


def test_P7_9_steady_state_of_a_random_walk_filter():
    q, r = 1e-4, 1e-2
    Pss = (-q + np.sqrt(q * q + 4 * q * r)) / 2
    Pk = 1.0
    for _ in range(2000):
        Pk = 1 / (1 / (Pk + q) + 1 / r)
    assert Pk == pytest.approx(Pss, rel=1e-10)


def test_zz_is_the_pauli_used_for_correlated_dephasing():
    assert np.allclose(X.pauli2(3, 3), np.kron(Z, Z))
