"""The hardware results regenerated from the shipped raw counts must equal the saved results of the research project
EXACTLY (relative tolerance 1e-11: the same arithmetic, only summation order may differ). This covers every Round 1 /
Round 2 fidelity with its bootstrap sigma, the readout calibration, the TOMO-9 / TOMO-12 matrices and learned
corrections, the detector, the per-circuit check, the comparison with the predictions, H1-H11, the twirl model, the
recovered fraction with its Fieller intervals and bootstrap (Part 9.8, P8.10), the error budget and the N2 analysis
(8e), and the job timing. Claim O8.12 of Part 9.8: every cached 8d number regenerated from the raw counts."""
import pytest

from telecorr.compare import assert_tree_close
from telecorr.hardware import analysis as A
from telecorr.hardware import data as D
from telecorr.hardware.__main__ import check


@pytest.fixture(scope="module")
def out():
    return A.reproduce()


@pytest.fixture(scope="module")
def ref():
    return {"R1": D.reference("8d_round1.json"), "R2": D.reference("8d_round2.json"), "CMP": D.reference("8d_compare.json")}


def test_round1_equals_saved(out, ref):
    n = assert_tree_close(out["r1"], ref["R1"]["r1"], min_leaves=1500)
    assert n >= 1500


def test_round2_equals_saved(out, ref):
    assert_tree_close(out["r2"], ref["R2"]["r2"], min_leaves=800)


def test_corrections_and_n0_rule(out, ref):
    assert_tree_close(out["corr"], ref["R1"]["corr"])
    assert out["corr"]["fingerprint_angle"] > 0.05 and out["corr"]["N0"] is not None   # N0 block was run
    # Round 2 was built with the corrections that the addendum recorded before sending
    assert_tree_close(out["corr"], D.preregistration("predictions_round2_addendum.json")["meta"]["corrections"])


def test_per_circuit_checks(out, ref):
    assert_tree_close(out["check_round1"], ref["R1"]["check"])
    assert_tree_close(out["check_round2"], ref["R2"]["check"])


def test_compare_and_hypotheses(out, ref):
    assert_tree_close(out["compare"], ref["CMP"]["rows"])
    assert_tree_close(out["hypotheses"], ref["CMP"]["hypotheses"])


def test_twirl_model_equals_preregistered_file(out):
    tm = D.preregistration("predictions_round2_tomo_model.json")
    assert_tree_close(out["twirl_model"]["blocks"], tm["blocks"])
    assert_tree_close(out["twirl_model"]["check_round1"], tm["check_round1"])


def test_all_saved_results_match(out):
    for name, n, nbad in check(out):
        assert nbad == 0, name
        assert n > 0, name


def test_headline_numbers(out):
    """The numbers quoted in the paper / README (rounded as printed there)."""
    rows = {r["state"]: r for r in out["recovered_fraction"]["rows"]}
    assert round(rows["product_HT"]["f"], 3) == 1.057 and round(rows["bell_30_70"]["f"], 3) == 0.975
    assert [round(x, 3) for x in rows["product_HT"]["fieller95"]] == [0.952, 1.177]
    assert [round(x, 3) for x in rows["bell_30_70"]["fieller95"]] == [0.934, 1.017]
    for st in ("product_HT", "bell_30_70"):
        lo, hi = rows[st]["fieller95"]
        assert lo <= 1 <= hi                                       # H2: the whole injected loss recovered
    t = out["r1"]["tomo"]
    assert round(t["R1_nativeA_product_HT"]["F"], 4) == 0.8751 and round(t["R1_nativeA_bell_30_70"]["F"], 4) == 0.8608
