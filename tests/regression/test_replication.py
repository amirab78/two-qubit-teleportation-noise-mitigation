"""The replication of 25 October 2026 (11 jobs on ibm_kingston and ibm_fez), regenerated from the shipped raw counts,
must equal the saved per-job analyses of the research project EXACTLY (relative tolerance 1e-11), and model 2 (the
Pauli-twirl model) recomputed from the counts must equal the values in the pre-registered prediction files."""
import math

import pytest

from telecorr.compare import assert_tree_close
from telecorr.hardware import data as D
from telecorr.hardware import replication as X
from telecorr.hardware.__main__ import check_replication


@pytest.fixture(scope="module")
def out():
    return X.reproduce()


def blocks(out, job):
    return {b["block"]: b for b in out[job]["blocks"]}


def test_all_saved_results_match(out):
    total = 0
    for name, n, nbad in check_replication(out):
        assert nbad == 0, name
        assert n > 0, name
        total += n
    assert total >= 4000


def test_round1_equals_saved(out):
    assert_tree_close(out["K1"]["r1"], D.replication_reference("x_K1_ibm_kingston.json")["r1"], min_leaves=600)
    assert_tree_close(out["F1"]["r1"], D.replication_reference("x_F1_ibm_fez.json")["r1"], min_leaves=280)


def test_round2_was_built_from_round1(out):
    """K2 / F2 ran the corrections that the prediction files recorded before sending, and those are the
    corrections Round 1's counts give."""
    assert_tree_close(out["K1"]["corr"], D.replication_preregistration("x_pred_K234_ibm_kingston.json")["meta"]["corrections"])
    assert_tree_close(out["F1"]["corr"], D.replication_preregistration("x_pred_F2_ibm_fez.json")["meta"]["corrections"])


def test_replication_on_kingston(out):
    """Three weeks later on the same device: the fresh correction recovers >= 90 % of N1's loss and reaches the
    device alone within 2 sigma; the correction stored from 4 October still works."""
    b = blocks(out, "K2"); t1 = out["K1"]["r1"]["tomo"]
    assert [round(b[f"R2_N1corr_{s}"]["F"], 3) for s in ("product_HT", "bell_30_70")] == [0.889, 0.856]
    assert [round(t1[f"R1_nativeA_{s}"]["F"], 3) for s in ("product_HT", "bell_30_70")] == [0.891, 0.857]
    assert [round(b[f"R2_N1corr_{s}"]["recovered_fraction"], 3) for s in ("product_HT", "bell_30_70")] == [0.988, 0.998]
    assert [round(t1[f"R1_N1old_{s}"]["F"], 3) for s in ("product_HT", "bell_30_70")] == [0.879, 0.865]
    assert all(c["ok"] for c in out["K2"]["criteria"] if c["criterion"].startswith("REPLICATION"))


def test_second_device_without_tuning(out):
    """ibm_fez, new chain, same pipeline: N1 0.706 / 0.511 -> 0.910 / 0.891; both pre-registered criteria pass."""
    b = blocks(out, "F2"); t1 = out["F1"]["r1"]["tomo"]
    assert [round(t1[f"R1_N1none_{s}"]["F"], 3) for s in ("product_HT", "bell_30_70")] == [0.706, 0.511]
    assert [round(b[f"R2_N1corr_{s}"]["F"], 3) for s in ("product_HT", "bell_30_70")] == [0.910, 0.891]
    assert all(c["ok"] for c in out["F2"]["criteria"])


def test_waits_linear_and_corrected(out):
    """The waiting phase grows linearly (axis -z, both legs near 15 kHz today); a fresh TOMO-9 recovers more than
    3 sigma at every wait, e.g. 0.261 -> 0.785 at 16 us."""
    w = out["wait_khz"]
    assert round(w["q2"], 1) == 15.8 and round(w["q3"], 1) == 15.1
    assert max(abs(r) for q in (2, 3) for r in w[f"residuals_q{q}"]) < 0.03       # rad, against angles up to 1.6
    b = blocks(out, "K2")
    assert [round(b[f"R2_N2corr{t}_product_HT"]["gain"], 3) for t in ("4us", "8us", "16us")] == [0.069, 0.225, 0.524]
    assert round(out["K1"]["r1"]["tomo"]["R1_N2none16us_product_HT"]["F"], 3) == 0.261
    assert round(b["R2_N2corr16us_product_HT"]["F"], 3) == 0.785
    assert all(c["ok"] for c in out["K2"]["criteria"] if c["criterion"].startswith("WAIT"))


def test_dd_spectator_difference_and_absolute_failure(out):
    """Part 9.8's spectator phase holds as a DIFFERENCE (K3 - K4 = pi zeta t) but not absolutely: without spectator
    pulses (K4) DD still leaves about +0.09 / +0.08 rad where the model predicts 0."""
    z3 = {r["leg"]: r for r in out["K3"]["rotation"]}; z4 = {r["leg"]: r for r in out["K4"]["rotation"]}
    diff = [z3[q]["z measured (rad)"] - z4[q]["z measured (rad)"] for q in ("q2", "q3")]
    pred = [z3[q]["z model 2 (rad)"] for q in ("q2", "q3")]
    assert [round(d, 3) for d in diff] == [0.108, 0.111] and [round(p, 3) for p in pred] == [0.114, 0.127]
    assert [round(z4[q]["z measured (rad)"], 2) for q in ("q2", "q3")] == [0.09, 0.08]
    assert all(z4[q]["pull"] > 4 for q in ("q2", "q3"))                       # the absolute prediction fails
    b3, b4 = blocks(out, "K3")["R3_DDcorr_product_HT"], blocks(out, "K4")["R4_DDUcorr_product_HT"]
    assert round(b4["F"] - b3["F"], 3) == 0.059 and round(math.hypot(b3["sigma"], b4["sigma"]), 3) == 0.011


def test_g1_same_job(out):
    """In one job on Fez: product_HT's 'corrected above native' was drift, bell_30_70's is real, and removing Fez's
    own rotation (fingerprint correction) adds about 0.04 to the Bell state."""
    t = {(r["state"], r["comparison"][:2]): r for r in out["G1"]["tests"]}
    assert round(t[("product_HT", "N1")]["D"], 3) == 0.002 and abs(t[("product_HT", "N1")]["z"]) < 1
    assert round(t[("bell_30_70", "N1")]["D"], 3) == 0.024 and t[("bell_30_70", "N1")]["z"] > 2.5
    assert round(t[("bell_30_70", "fi")]["D"], 3) == 0.039 and t[("bell_30_70", "fi")]["z"] > 4


def test_g2_dd_rotation_grows_with_the_wait(out):
    """z(t) = a + b t under DD: the slope is about 8 sigma from zero on both legs (+2.5 / +2.7 kHz), so CONST is
    ruled out; q2 reads LINEAR, q3 is undecided (intercept -0.065 +/- 0.020)."""
    f = out["G2"]["fits"]
    khz = [f[q]["b_rad_per_us"] / (2 * math.pi) * 1e3 for q in ("q2", "q3")]
    assert [round(k, 1) for k in khz] == [2.5, 2.7]
    assert all(f[q]["b_rad_per_us"] / f[q]["b_sigma"] > 8 for q in ("q2", "q3"))
    assert not f["q2"]["CONST favoured (b ~ 0)"] and f["q2"]["LINEAR favoured (a ~ 0)"]
    assert not f["q3"]["CONST favoured (b ~ 0)"] and not f["q3"]["LINEAR favoured (a ~ 0)"]


def test_g3_correction_under_dd_barely_helps(out):
    """G3-A fails at both waits; DD + correction is worse than K2's fresh correction without DD by > 4 sigma."""
    c = {x["criterion"][:9]: x for x in out["G3"]["criteria"]}
    assert not c["G3-A 8us:"]["ok"] and not c["G3-A 16us"]["ok"]
    assert [round(c[k]["value"], 3) for k in ("G3-A 8us:", "G3-A 16us")] == [-0.002, 0.025]
    assert [round(c[k]["z"], 1) for k in ("G3-B 8us:", "G3-B 16us")] == [-6.6, -4.9]


def test_jobs_within_their_caps(out):
    rows = out["jobs"]
    assert len(rows) == 11 and sum(r["usage_s"] for r in rows) == 318
    assert all(r["usage_s"] <= r["cap_s"] for r in rows)
    assert all(abs(r["usage_s"] - r["estimate_detailed_s"]) < 3 for r in rows)
