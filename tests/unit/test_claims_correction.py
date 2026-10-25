"""Proven results of steps 3-5 (Parts 9.3-9.5): probe cost, its landscape, the ceiling and Rotosolve."""
import numpy as np
import pytest
from scipy.optimize import minimize

from telecorr import calibrate as K
from telecorr import channels as C
from telecorr import probes as Pr
from telecorr import protocol as P
from telecorr.ceiling import ceiling, procrustes
from telecorr.paulis import I2
from telecorr.rotations import angles_of, axis_rotation, correction_rotations, so3, u_matrix
from telecorr.states import message_state
from telecorr.twirl import f_avg_of

RNG = np.random.default_rng(93)


def random_rotation(rng):
    return so3(u_matrix(*rng.uniform(-np.pi, np.pi, 3)))


def random_T(rng):
    fam = rng.choice(list(C.FAMILIES))
    D = C.incoherent_T(**C.leg_parameters(fam, rng.uniform(0, 0.3), 2))
    return D @ random_rotation(rng)


def test_T3_5_probe_probabilities_are_products_of_leg_factors():
    """All four outcome probabilities of every probe circuit: twirl formula = the qiskit circuit (Aer, density matrix)."""
    for _ in range(2):
        x = RNG.uniform(-np.pi, np.pi, 6)
        legs = C.leg_kraus(C.random_noise(RNG), "M", 0.1)
        T = {q: C.T_matrix(legs[q]) for q in (2, 3)}
        for name in Pr.PROBE_AXIS:
            rho = P.aer_bob_state(P.circuit(name, "unitary", unprepare=True, correction=x), C.pair_kraus(legs))
            assert np.allclose(np.real(np.diag(rho)), Pr.probe_probs(x, T, (name,))[0], atol=1e-12)


def test_T3_6_T5_9_the_cost_is_sandwiched_by_the_average_fidelity():
    for _ in range(300):
        T = {2: random_T(RNG), 3: random_T(RNG)}
        x = RNG.uniform(-np.pi, np.pi, 6)
        c = Pr.cost(x, T)
        lo, hi = Pr.sandwich(f_avg_of(correction_rotations(x), T))
        assert lo - 1e-12 <= c <= hi + 1e-12


def test_T3_8_one_probe_is_blind_to_rotations_about_its_own_axis():
    T = {2: so3(axis_rotation(0, 1.0)), 3: so3(axis_rotation(0, 2.0))}       # residual rotations about x
    assert Pr.cost(None, T, probes=("|++>",)) == pytest.approx(0, abs=1e-12)
    assert f_avg_of(correction_rotations(None), T) < 0.6


def test_P4_1_cost_is_bilinear_in_the_squared_quaternions():
    for _ in range(50):
        q = [RNG.normal(size=4) for _ in range(2)]
        q = [v / np.linalg.norm(v) for v in q]
        W = [v[0] * I2 - 1j * sum(c * s for c, s in zip(v[1:], C.SIGMA)) for v in q]
        T = {2: so3(W[0]), 3: so3(W[1])}
        assert Pr.cost(None, T) == pytest.approx(Pr.quaternion_cost(q[0] ** 2, q[1] ** 2), abs=1e-12)


def test_T4_2_no_spurious_local_minimum():
    """Along the path u(t) = (1 - t) u + t e0 (both legs) the cost decreases strictly to 0 from every point."""
    e0 = np.array([1.0, 0, 0, 0]); t = np.linspace(0, 1, 41)
    for _ in range(500):
        u, v = RNG.dirichlet(np.ones(4) * 0.5), RNG.dirichlet(np.ones(4) * 0.5)
        c = [Pr.quaternion_cost((1 - s) * u + s * e0, (1 - s) * v + s * e0) for s in t]
        assert np.all(np.diff(c) < 1e-15) and c[-1] == pytest.approx(0, abs=1e-15)


def test_T4_3_critical_values():
    e = np.eye(4)
    cases = {0: (e[0], e[0]), 2 / 3: (e[0], e[1]), 5 / 6: (e[1], (e[1] + e[2]) / 2), 8 / 9: (e[1], (e[1] + e[2] + e[3]) / 3),
             1: (e[1], e[2])}
    for value, (u, v) in cases.items():
        assert Pr.quaternion_cost(u, v) == pytest.approx(value, abs=1e-15)


def test_P4_10_an_accepted_start_guarantees_a_high_fidelity():
    """An exact cost C < 0.01 implies F_e > 0.985, a lower bound on the fidelity of every message."""
    for _ in range(2000):
        T = {q: so3(u_matrix(*RNG.normal(0, 0.05, 3))) for q in (2, 3)}
        c = Pr.cost(None, T)
        pI = {q: (1 + np.trace(T[q])) / 4 for q in (2, 3)}
        assert 1 - pI[2] * pI[3] <= 1.5 * c + 1e-12
    assert 1 - 1.5 * 0.01 == pytest.approx(0.985)


def test_T5_4_average_fidelity_formula():
    for _ in range(4):
        legs = C.leg_kraus(C.random_noise(RNG), RNG.choice(list(C.FAMILIES)), RNG.uniform(0, 0.3))
        T = {q: C.T_matrix(legs[q]) for q in (2, 3)}
        assert P.average_fidelity(legs=legs) == pytest.approx(f_avg_of(correction_rotations(None), T), abs=1e-12)


def test_T5_6_ceiling_is_procrustes_and_the_inverse_for_symmetric_noise():
    for _ in range(100):
        T = RNG.normal(size=(3, 3))
        R = procrustes(T)
        s = np.linalg.svd(T, compute_uv=False)
        assert np.allclose(R.T @ R, np.eye(3)) and np.linalg.det(R) == pytest.approx(1)
        assert np.trace(R @ T) == pytest.approx(s[0] + s[1] + np.sign(np.linalg.det(T)) * s[2], abs=1e-12)
        assert max(np.trace(random_rotation(RNG) @ T) for _ in range(50)) <= np.trace(R @ T) + 1e-12
    for fam in C.FAMILIES:                                   # Stage 5's channels: the ceiling is the inverse rotation
        T = C.leg_T(None, fam, 0.3)
        Rinv = {q: so3(C.coherent_unitaries()[q]).T for q in (2, 3)}
        assert ceiling(T)[1] == pytest.approx(f_avg_of(Rinv, T), abs=1e-12)


def test_P5_10_probe_bias_is_real_but_second_order():
    """Amplitude damping p = 0.3: the 3-probe optimum misses the ceiling by about 1e-6 (values of the study)."""
    T = C.leg_T(None, "A", 0.3)
    x0, Fc = ceiling(T)
    best = min((minimize(lambda x: Pr.cost(x, T), x0 + RNG.normal(0, 0.05, 6), method="BFGS", options={"gtol": 1e-12})
                for _ in range(4)), key=lambda r: r.fun)
    deficit = Fc - f_avg_of(correction_rotations(best.x), T)
    R = correction_rotations(best.x)
    resid = max(np.arccos(np.clip((np.trace(R[q] @ so3(C.coherent_unitaries()[q])) - 1) / 2, -1, 1)) for q in (2, 3))
    assert deficit == pytest.approx(1.194e-6, rel=0.02)
    assert resid == pytest.approx(2.486e-3, rel=0.02)


def test_P5_15_rotosolve_coordinates_are_exact_sinusoids():
    T = {2: random_T(RNG), 3: random_T(RNG)}
    x = RNG.uniform(-np.pi, np.pi, 6)
    for d in range(6):
        f = lambda t: Pr.cost(K.recentre(x, np.eye(6)[d] * t), T)            # noqa: E731
        c0, cp, cm = f(0), f(np.pi / 2), f(-np.pi / 2)
        a, b, s = (cp + cm) / 2, c0 - (cp + cm) / 2, (cp - cm) / 2
        for t in RNG.uniform(-np.pi, np.pi, 5):
            assert f(t) == pytest.approx(a + b * np.cos(t) + s * np.sin(t), abs=1e-12)
        t_star = np.arctan2(cm - cp, cp + cm - 2 * c0)
        grid = np.linspace(-np.pi, np.pi, 721)
        assert f(t_star) <= min(f(t) for t in grid) + 1e-12


def test_rotosolve_with_exact_costs_descends_monotonically():
    class Exact(K.ProbeSampler):
        def cost(self, params, shots):
            return Pr.cost(params, self.T)
    T = {2: random_T(RNG), 3: random_T(RNG)}
    s = Exact(T, RNG); trace = []
    K.rotosolve(np.zeros(6), s, phases=((1000, 3, 3),), trace=trace)
    costs = [Pr.cost(x, T) for _, x in trace]
    assert np.all(np.diff(costs) < 1e-12)


def test_recentre_keeps_the_correction_unchanged_at_zero():
    x = RNG.uniform(-np.pi, np.pi, 6)
    y = K.recentre(x, np.zeros(6))
    for i in (0, 3):
        V, W = u_matrix(*x[i:i + 3]), u_matrix(*y[i:i + 3])
        assert abs(abs(np.trace(V.conj().T @ W)) - 2) < 1e-12
    assert np.allclose(angles_of(u_matrix(*y[:3])), y[:3])
