"""Command line: rebuild every hardware result from the shipped raw counts.

Roadmap: step 8 (real hardware).

    python -m telecorr.hardware reproduce [--out results/hardware_reproduced.json] [--check]

``--check`` compares the result with the saved results of the research project (``data/reference``) and exits
with status 1 if any value differs.
"""
from __future__ import annotations

import argparse
import json
import sys
import time

from ..compare import tree_diff
from . import analysis as A
from . import data as D


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


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="python -m telecorr.hardware", description="hardware analysis from the raw counts")
    ap.add_argument("command", choices=["reproduce"])
    ap.add_argument("--out", default=None, help="write the full result as JSON")
    ap.add_argument("--check", action="store_true", help="compare with the saved results")
    a = ap.parse_args(argv)
    t0 = time.time()
    out = A.reproduce()
    print(f"hardware analysis done in {time.time() - t0:.1f} s")
    if a.out:
        with open(a.out, "w", encoding="utf8") as fh:
            json.dump(out, fh, indent=1, default=float)
        print("written", a.out)
    if a.check:
        ok = True
        for name, n, nbad in check(out):
            print(f"  {'OK ' if nbad == 0 else 'DIFF'} {name}: {n} values, {nbad} differ")
            ok &= nbad == 0
        return 0 if ok else 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
