"""Saved simulation results of steps 1-7 reproduced by the package.

Deterministic results are compared exactly (to 1e-9 or to the printed digits). Results of a shot-noise calibration
cannot match a saved run draw by draw (the study sampled with Aer and its own random streams), so those tests check
the stated property with a tolerance written next to it. The Stage 7 tracking study used numpy streams only and is
reproduced exactly.
"""
import numpy as np
import pytest

from telecorr import calibrate as K
from telecorr import channels as C
from telecorr import correlated as X
from telecorr import probes as Pr
from telecorr import protocol as P
from telecorr import tracking as Tr
from telecorr.ceiling import ceiling
from telecorr.reference import number, simulation
from telecorr.states import MESSAGES, message_state

MSGS = [message_state(n) for n in MESSAGES]


def test_step1_unpreparation_reads_00_in_every_shot():
    for name in MESSAGES:
        assert P.bob_counts(P.circuit(name, "measured", unprepare=True), 8192, seed=1234) == {"00": 8192}


def test_step2_fidelity_and_purity_of_the_coherent_channel():
    U = C.coherent_unitaries()
    for name, kF, kp in (("product_HT", "n0003", "n0005"), ("bell_30_70", "n0004", "n0006")):
        rho = P.teleport(message_state(name), legs={2: [U[2]], 3: [U[3]]})["bob"]
        assert P.fidelity(message_state(name), rho) == pytest.approx(number(kF), abs=1e-12)
        assert P.purity(rho) == pytest.approx(number(kp), abs=1e-12)


def test_step2_best_bob_side_correction_reaches_lambda_max():
    U = C.coherent_unitaries()
    rho = P.teleport(message_state("product_HT"), legs={2: [U[2]], 3: [U[3]]})["bob"]
    assert np.linalg.eigvalsh(rho)[-1] == pytest.approx(0.7963045758580273, abs=1e-12)


def test_step3_saved_stage3_correction():
    """Stage 3's saved angles: probe cost 8.93e-5, F = 0.999906 / 0.999874, residual 0.0164 / 0.0163 rad."""
    x = np.array(simulation("stage3_calibration.json")["params"])
    T = C.leg_T()
    s = K.scores(x, T, MSGS)
    assert [round(f, 6) for f in s["F"]] == [float(number("n0007")), float(number("n0008"))]
    assert [round(a, 4) for a in K.residual_angles(x, T)] == [float(number("n0103")), float(number("n0104"))]
    assert Pr.cost(x, T) == pytest.approx(8.934832523888048e-05, rel=1e-9)


def test_step3_fresh_calibration_reaches_the_shot_noise_floor():
    """A fresh COBYLA calibration from counts (seed fixed): F_avg > 0.995 (30 seeds: median 1 - F_avg = 4.6e-4)."""
    T = C.leg_T()
    x = K.cobyla(np.zeros(6), K.ProbeSampler(T, np.random.default_rng(2024)))
    assert K.scores(x, T, MSGS)["F_avg"] > 0.995


def test_step4_hard_case():
    H = C.leg_T(C.hard_case())
    assert [round(f, 6) for f in K.scores(None, H, MSGS)["F"]] == [float(number("n0011")), float(number("n0012"))]
    auto = K.AutoCalibrator(K.ProbeSampler(H, np.random.default_rng(5)), seed=5)
    auto.step()
    assert min(K.scores(auto.params, H, MSGS)["F"]) > 0.985          # the accept-rule guarantee (P4.10)


def test_step5_realistic_channel_E4():
    T = C.leg_T(None, "M", 0.1)
    assert K.scores(None, T, MSGS)["F"] == pytest.approx([number("n0019"), number("n0020")], abs=1e-12)
    xc, _ = ceiling(T)
    assert [round(f, 6) for f in K.scores(xc, T, MSGS)["F"]] == [float(number("n0023")), float(number("n0024"))]
    cal = K.RelativeAutoCalibrator(K.ProbeSampler(T, np.random.default_rng(5)), seed=5)
    cal.step()
    gap = max(c - f for c, f in zip(K.scores(xc, T, MSGS)["F"], K.scores(cal.params, T, MSGS)["F"]))
    assert gap < 1e-4                                                # saved run: 4e-6 / 3e-6


def test_step6_nominal_correlated_channel_table():
    saved = {r["label"]: r for r in simulation("6a_nominal.json")["rows"]}
    T = X.T15()
    W = X.ceilings(T)
    for mine, theirs in (("no correction", "no correction"), ("per-leg Procrustes (TOMO-9)", "per-leg Procrustes (= TOMO-9 on each leg)"),
                         ("best local", "best local (6 parameters)"), ("best two-qubit", "best two-qubit (15)")):
        r = saved[theirs]
        assert X.favg(W[mine], T) == pytest.approx(r["favg"], abs=1e-9)
        assert X.fidelity(MSGS[0], W[mine], T) == pytest.approx(r["product_HT"], abs=1e-9)
        assert X.fidelity(MSGS[1], W[mine], T) == pytest.approx(r["bell"], abs=1e-9)
    assert X.favg_bound(T) == pytest.approx(simulation("6a_nominal.json")["bound"], abs=1e-12)


@pytest.mark.slow
def test_step7_tracking_slow_drift_exactly():
    """Stage 7f, scenario 'slow' (200 episodes, its seeds): ABC and TOMO-9 every step reproduce the saved values."""
    saved = simulation("7f_tracking.json")["slow"]
    Ts = Tr.drift_scenarios(np.random.default_rng(102), 200, **Tr.SCENARIOS["slow"])
    for tg in (3e-4, 1e-4):
        r = Tr.summary(Tr.run_abc(Ts, tg, seed=111))
        s = saved[f"ABC target {tg:g}"]
        assert r["mean_gap"] == pytest.approx(s["mean_gap"], rel=1e-9)
        assert r["shots_per_step"] == pytest.approx(s["shots_per_step"], rel=1e-12)
    r = Tr.summary(Tr.run_tomo9(Ts, 1000, seed=121))
    assert r["mean_gap"] == pytest.approx(saved["TOMO-9 every step, N = 1000"]["mean_gap"], rel=1e-9)
