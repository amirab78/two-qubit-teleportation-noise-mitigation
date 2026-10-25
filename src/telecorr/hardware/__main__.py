"""Command line: rebuild every hardware result from the shipped raw counts.

Roadmap: step 8 (real hardware).

    python -m telecorr.hardware reproduce [--run all|2026-10-04|2026-10-25] [--out results/hardware_reproduced.json] [--check]

``--run`` picks the experiment: the 4 October jobs, the replication of 25 October, or both (default).
``--check`` compares the result with the saved results of the research project (``data/reference``) and the
pre-registered predictions, and exits with status 1 if any value differs.
"""
from __future__ import annotations

import argparse
import json
import sys
import time

from ..compare import tree_diff
from . import analysis as A
from . import data as D
from . import replication as X


def check(out: dict) -> list[tuple[str, int, int]]:
    """(name, values compared, mismatches) for every saved hardware result."""
    R1, R2 = D.reference("8d_round1.json"), D.reference("8d_round2.json")
    cmp = D.reference("8d_compare.json")
    tm = D.preregistration("predictions_round2_tomo_model.json")
    rec = next(c for c in D.reference("9_8_checks.json")["claims"] if c["key"] == "recov")
    rows98 = next(ch for ch in rec["checks"] if "rows" in ch)["rows"]
    boot = D.reference("9_8_boot.json")["results"]
    dd98 = next(ch for ch in next(c for c in D.reference("9_8_checks.json")["claims"] if c["key"] == "ddobs")["checks"]
                if "rows" in ch)
    keys = ("DD", "TOMO+DD", "chi2", "p")
    dd_mine = {r["scenario"][:3]: {k: r[k] for k in keys} for r in out["dd"]["rows"]}
    dd_ref = {r["scenario"][:3]: {k: r[k] for k in keys} for r in dd98["rows"]}
    pairs = [("Round 1 (8d_round1.json: r1)", out["r1"], R1["r1"]),
             ("Round 2 (8d_round2.json: r2)", out["r2"], R2["r2"]),
             ("Round 2 corrections", out["corr"], R1["corr"]),
             ("per-circuit check, Round 1", out["check_round1"], R1["check"]),
             ("per-circuit check, Round 2", out["check_round2"], R2["check"]),
             ("hardware vs predictions (8d_compare.json)", out["compare"], cmp["rows"]),
             ("hypotheses H1-H11", out["hypotheses"], cmp["hypotheses"]),
             ("twirl model (predictions_round2_tomo_model.json)", out["twirl_model"]["blocks"], tm["blocks"]),
             ("recovered fraction, Fieller (Part 9.8)", [{k: r[k] for k in rows98[0] if k in r} for r in out["recovered_fraction"]["rows"]],
              [{k: r[k] for k in r if k in out["recovered_fraction"]["rows"][0]} for r in rows98]),
             ("recovered fraction, bootstrap (9_8_boot.json)", out["recovered_fraction"]["bootstrap"]["results"],
              {st: {k: v[k] for k in ("f_q", "f_sd", "P_f_gt_1", "first10")} for st, v in boot.items()}),
             ("error budget (8e_error_budget.json)", out["error_budget"], D.reference("8e_error_budget.json")),
             ("waiting channel N2 (8e_n2_analysis.json)", out["n2"], D.reference("8e_n2_analysis.json")),
             ("DD spectator scenarios (Part 9.8, O8.7)", [dd_mine, out["dd"]["spectator_rad"]], [dd_ref, dd98["spectator_rad"]]),
             ("job timing (8d_jobs_table.json)", out["jobs"],
              [{k: r[k] for k in out["jobs"][0] if k in r} for r in D.reference("8d_jobs_table.json")])]
    res = []
    for name, mine, saved in pairs:
        n, bad = tree_diff(mine, saved)
        res.append((name, n, len(bad)))
    return res


REPLICATION_FILES = {"K0": "x_K0_ibm_kingston.json", "K1": "x_K1_ibm_kingston.json", "K2": "x_K2_ibm_kingston.json",
                     "K3": "x_K3_ibm_kingston.json", "K4": "x_K4_ibm_kingston.json", "F0": "x_F0_ibm_fez.json",
                     "F1": "x_F1_ibm_fez.json", "F2": "x_F2_ibm_fez.json", "G1": "x_G1_ibm_fez.json",
                     "G2": "x_G2_ibm_kingston.json", "G3": "x_G3_ibm_kingston.json"}
TIMING = ("usage_s", "estimate_detailed_s", "cap_s", "queue_s", "run_wall_s")


def check_replication(out: dict) -> list[tuple[str, int, int]]:
    """(name, values compared, mismatches) for every saved result and prediction of the replication."""
    P = D.replication_preregistration
    k234, f2, g, g3 = (P("x_pred_K234_ibm_kingston.json"), P("x_pred_F2_ibm_fez.json"), P("x_pred_G.json"),
                       P("x_pred_G3.json"))
    pairs = []
    for job, name in REPLICATION_FILES.items():
        saved = D.replication_reference(name)
        keys = [k for k in out[job] if k in saved]
        pairs.append((f"{job} ({name})", {k: out[job][k] for k in keys}, {k: saved[k] for k in keys}))
    pairs += [("model 2, K2 (x_pred_K234_ibm_kingston.json)", out["K2"]["model_2"], k234["model_2"]["blocks"]),
              ("model 2, DD of K3 / K4 (x_pred_K234_ibm_kingston.json)", {"dd": out["K3"]["model_2"], "ddu": out["K4"]["model_2"]},
               k234["model_2"]["dd"]),
              ("model 2, F2 (x_pred_F2_ibm_fez.json)", out["F2"]["model_2"], f2["model_2"]["blocks"]),
              ("model 2 and gap hypotheses, G1 (x_pred_G.json)", [out["G1"]["model_2"], out["G1"]["gap_hypotheses"]],
               [g["G1"]["model_2"], g["G1"]["gap_hypotheses"]]),
              ("model 2, corrections and references, G3 (x_pred_G3.json)",
               [out["G3"]["model_2"], out["G3"]["references"], out["G2"]["corrections_for_G3"]],
               [g3["model_2"], g3["references"], g3["meta"]["corrections_from_G2"]]),
              ("K1 corrections and waiting frequency (x_pred_K234_ibm_kingston.json, x_K1_retro_wait.json)",
               [out["K1"]["corr"], out["wait_khz"], out["wait_khz"]],
               [k234["meta"]["corrections"], k234["meta"]["wait_detuning_K1_kHz"],
                D.replication_reference("x_K1_retro_wait.json")["wait_detuning_K1_kHz"]]),
              ("job timing (usage of every x_<job>.json)", [{k: r[k] for k in TIMING} for r in out["jobs"]],
               [{k: D.replication_reference(REPLICATION_FILES[r["job"]])["usage"][k] for k in TIMING} for r in out["jobs"]])]
    res = []
    for name, mine, saved in pairs:
        n, bad = tree_diff(mine, saved)
        res.append((name, n, len(bad)))
    return res


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="python -m telecorr.hardware", description="hardware analysis from the raw counts")
    ap.add_argument("command", choices=["reproduce"])
    ap.add_argument("--run", choices=["all", "2026-10-04", "2026-10-25"], default="all",
                    help="the 4 October experiment, the replication of 25 October, or both")
    ap.add_argument("--out", default=None, help="write the full result as JSON")
    ap.add_argument("--check", action="store_true", help="compare with the saved results")
    a = ap.parse_args(argv)
    t0 = time.time()
    out = {}
    if a.run in ("all", "2026-10-04"):
        out["2026-10-04"] = A.reproduce()
    if a.run in ("all", "2026-10-25"):
        out["2026-10-25"] = X.reproduce()
    print(f"hardware analysis done in {time.time() - t0:.1f} s")
    if a.out:
        with open(a.out, "w", encoding="utf8") as fh:
            json.dump(out, fh, indent=1, default=float)
        print("written", a.out)
    if a.check:
        ok, total = True, 0
        for run, fn in (("2026-10-04", check), ("2026-10-25", check_replication)):
            if run not in out:
                continue
            print(f"{run}:")
            for name, n, nbad in fn(out[run]):
                print(f"  {'OK ' if nbad == 0 else 'DIFF'} {name}: {n} values, {nbad} differ")
                ok &= nbad == 0; total += n
        print(f"{total} values compared, {'none differ' if ok else 'SOME DIFFER'}")
        return 0 if ok else 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
