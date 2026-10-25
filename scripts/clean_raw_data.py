"""Build ``data/`` from the original research project (run once by the authors; kept to document how the
shipped data were made).

usage:  python scripts/clean_raw_data.py --project PATH_TO_RESEARCH_PROJECT

What it does
------------
* hardware counts of the 4 IBM jobs -> ``counts/<job id>.json.gz``: the counts and per-shot outcomes are kept
  unchanged; from the job record and metrics only fields that describe the experiment are kept (no account,
  instance, quota or billing fields; the usage-detail metric id is dropped);
* the calibration record in force when the jobs ran (identical for the 4 jobs) -> one gzip file, with every
  time stamp converted from the local time zone of the recording computer to UTC;
* frozen circuits (QPY + parameters + label index), predictions, their SHA-256 records and the approvals are
  copied byte for byte, so the recorded hashes can be checked;
* saved results used as regression targets -> ``data/reference/`` (simulation results in ``simulation/``) (byte for byte, except that account quota
  fields and the usage-detail metric id are removed from the two 8d round files);
* ``data/MANIFEST.json``: every shipped file with its source path, source SHA-256, shipped SHA-256 and the
  transformation applied; then a check that the cleaned counts equal the original counts, pub by pub.
"""
from __future__ import annotations

import argparse
import datetime as dt
import gzip
import hashlib
import io
import json
import os
import re
import shutil

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
DATA = os.path.join(REPO, "data")
HW = os.path.join(DATA, "hardware", "ibm_kingston_2026-10-04")

JOBS = {  # job id -> (role, frozen job name)
    "db16kn2vog1s73fhnpog": ("pilot (Round 0)", "R0_main_ibm_kingston"),
    "db16p1bid5ic73eqnsug": ("Round 1", "R1_main_ibm_kingston"),
    "db175bjid5ic73eqoav0": ("Round 2 main", "R2_main_ibm_kingston"),
    "db175givog1s73fhocjg": ("Round 2 dynamical decoupling", "R2_dd_ibm_kingston"),
}
FROZEN_DIR = {"R0": "Stage8_Hardware/8c_Twin_Study/results/jobs/ibm_kingston",
              "R1": "Stage8_Hardware/8c_Twin_Study/results/jobs/ibm_kingston",
              "R2": "Stage8_Hardware/8d_Hardware_Runs/results/jobs/ibm_kingston"}
CALIBRATION_STAMP = "20261004T1434Z"

RECORD_KEEP = ("id", "where", "backend", "round", "purpose", "labels", "n_pubs", "shots", "shots_total", "estimate_s",
               "max_seconds", "dd", "twirl", "sampler", "submitted_utc", "usage_s", "hash", "ibm_estimates_queued")
PREREG = [  # (source, shipped name)
    ("Stage8_Hardware/8c_Twin_Study/results/predictions_round1.json", "predictions_round1.json"),
    ("Stage8_Hardware/8c_Twin_Study/results/predictions_round2.json", "predictions_round2.json"),
    ("Stage8_Hardware/8c_Twin_Study/results/predictions_round1_ibm_fez.json", "predictions_round1_ibm_fez.json"),
    ("Stage8_Hardware/8c_Twin_Study/results/predictions_round2_ibm_fez.json", "predictions_round2_ibm_fez.json"),
    ("Stage8_Hardware/8c_Twin_Study/results/8c_predictions_sha256.json", "8c_predictions_sha256.json"),
    ("Stage8_Hardware/8d_Hardware_Runs/results/predictions_round2_addendum.json", "predictions_round2_addendum.json"),
    ("Stage8_Hardware/8d_Hardware_Runs/results/predictions_round2_tomo_model.json", "predictions_round2_tomo_model.json"),
    ("Stage8_Hardware/8d_Hardware_Runs/results/8d_predictions_sha256.json", "8d_predictions_sha256.json"),
    ("Stage8_Hardware/8c_Twin_Study/results/8c_choices.json", "8c_choices.json"),
    ("Stage8_Hardware/shared/jobs/approved.json", "approvals.json"),
]
REFERENCE = [  # (source, shipped name, cleaned?)
    ("Stage8_Hardware/8d_Hardware_Runs/results/8d_round1.json", "8d_round1.json", True),
    ("Stage8_Hardware/8d_Hardware_Runs/results/8d_round2.json", "8d_round2.json", True),
    ("Stage8_Hardware/8d_Hardware_Runs/results/8d_compare.json", "8d_compare.json", False),
    ("Stage8_Hardware/8d_Hardware_Runs/results/8d_jobs_table.json", "8d_jobs_table.json", False),
    ("Stage8_Hardware/8c_Twin_Study/results/8c_ladder.json", "8c_ladder.json", False),
    ("Stage8_Hardware/8d_Hardware_Runs/results/8d_postdiction.json", "8d_postdiction.json", False),
    ("Stage8_Hardware/8e_Results_Summary/results/8e_numbers.json", "8e_numbers.json", False),
    ("Stage8_Hardware/8e_Results_Summary/results/8e_error_budget.json", "8e_error_budget.json", False),
    ("Stage8_Hardware/8e_Results_Summary/results/8e_n2_analysis.json", "8e_n2_analysis.json", False),
    ("Stage8_Hardware/8e_Results_Summary/results/8e_twirl_check.json", "8e_twirl_check.json", False),
    ("Stage9_Theory_Reports/9_8_Stage8_Hardware/results/9_8_checks.json", "9_8_checks.json", False),
    ("Stage9_Theory_Reports/9_8_Stage8_Hardware/results/9_8_boot.json", "9_8_boot.json", False),
    ("Stage9_Theory_Reports/9_8_Stage8_Hardware/results/9_8_numbers.json", "9_8_numbers.json", False),
    ("Final_Report/results/final_numbers.json", "final_numbers.json", "numbers"),
]
SIMULATION = [  # saved simulation results used by notebooks 3, 6 and 7 (source, shipped name), byte copies
    ("Stage3_Correction/calibration.json", "stage3_calibration.json"),
    ("Stage6_Correlated_Noise/results/6a_nominal.json", "6a_nominal.json"),
    ("Stage6_Correlated_Noise/results/6a_tomo75_design.json", "6a_tomo75_design.json"),
    ("Stage6_Correlated_Noise/results/6b_detector_counts.json", "6b_detector_counts.json"),
    ("Stage7_Reinforcement_Learning/results/7f_tracking.json", "7f_tracking.json"),
    ("Stage7_Reinforcement_Learning/results/7g_controller.json", "7g_controller.json"),
]
QUOTA_KEYS = {"budget"}                     # account quota of the IBM plan: not part of the experiment
QUOTA_PATHS = ("usage_limit_seconds", "usage_remaining_seconds")


def sha256_bytes(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def sha256_file(path: str) -> str:
    with open(path, "rb") as fh:
        return sha256_bytes(fh.read())


def write_gz_json(path: str, obj) -> bytes:
    """Deterministic gzip of compact JSON (fixed mtime, no file name): same input -> same bytes."""
    raw = json.dumps(obj, separators=(",", ":")).encode("utf8")
    buf = io.BytesIO()
    with gzip.GzipFile(filename="", mode="wb", fileobj=buf, compresslevel=9, mtime=0) as gz:
        gz.write(raw)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "wb") as fh:
        fh.write(buf.getvalue())
    return buf.getvalue()


def clean_metrics(m: dict | None) -> dict | None:
    if not m:
        return m
    out = {k: m[k] for k in ("timestamps", "circuits_execution_time_ns", "qiskit_version") if k in m}
    if "usage" in m:
        out["usage"] = {k: v for k, v in m["usage"].items() if k != "details"}
    return out


def clean_job(raw: dict) -> dict:
    return {"id": raw["id"], "record": {k: raw["record"][k] for k in RECORD_KEEP if k in raw["record"]},
            "pubs": raw["pubs"], "metrics": clean_metrics(raw.get("metrics"))}


_TS = re.compile(r"^(\d{4}-\d\d-\d\dT\d\d:\d\d:\d\d(?:\.\d+)?)([+-]\d\d:\d\d)$")


def to_utc(obj):
    """Convert every ISO time stamp with an offset to UTC ('...Z'); the instant is unchanged."""
    if isinstance(obj, dict):
        return {k: to_utc(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [to_utc(v) for v in obj]
    if isinstance(obj, str) and _TS.match(obj):
        t = dt.datetime.fromisoformat(obj).astimezone(dt.timezone.utc)
        return t.isoformat().replace("+00:00", "Z")
    return obj


def clean_reference(obj):
    """Drop the quota fields and the usage-detail metric id from a saved 8d round file."""
    if isinstance(obj, dict):
        out = {}
        for k, v in obj.items():
            if k in QUOTA_KEYS:
                continue
            if k == "metrics_usage" and isinstance(v, dict):
                v = {a: b for a, b in v.items() if a != "details"}
            out[k] = clean_reference(v)
        return out
    if isinstance(obj, list):
        return [clean_reference(v) for v in obj]
    return obj


def clean_numbers(obj: dict) -> dict:
    """Drop the plan's quota (limit, remaining) from the final report's number log; QPU seconds used stay."""
    nums = {k: v for k, v in obj["numbers"].items() if not str(v.get("path", "")).endswith(QUOTA_PATHS)}
    return dict(obj, numbers=nums, n=len(nums))


def main(project: str) -> None:
    src = lambda rel: os.path.join(project, *rel.split("/"))            # noqa: E731
    rel = lambda p: os.path.relpath(p, REPO).replace(os.sep, "/")        # noqa: E731
    manifest = []
    if os.path.isdir(HW):
        shutil.rmtree(HW)
    ref_dir = os.path.join(DATA, "reference")
    if os.path.isdir(ref_dir):
        shutil.rmtree(ref_dir)

    # 1. counts of the 4 hardware jobs
    for jid, (role, name) in JOBS.items():
        s = src(f"Stage8_Hardware/shared/jobs/raw/{jid}.json")
        with open(s, "rb") as fh:
            b = fh.read()
        raw = json.loads(b)
        cleaned = clean_job(raw)
        assert cleaned["pubs"] == raw["pubs"], "counts must be unchanged"
        dst = os.path.join(HW, "counts", f"{jid}.json.gz")
        out = write_gz_json(dst, cleaned)
        manifest.append({"path": rel(dst), "what": f"counts and per-shot outcomes, {role} ({name})",
                         "source": f"Stage8_Hardware/shared/jobs/raw/{jid}.json", "source_sha256": sha256_bytes(b),
                         "sha256": sha256_bytes(out), "bytes": len(out),
                         "transform": "kept counts unchanged; record and metrics reduced to experiment fields; gzip"})

    # 2. the job-time calibration record (identical for the 4 jobs)
    props = [src(f"Stage8_Hardware/shared/snapshots/ibm_kingston/jobs/props_{j}.json") for j in JOBS]
    hashes = {sha256_file(p) for p in props}
    assert len(hashes) == 1, "the 4 job-time calibration records must be identical"
    with open(props[0], "rb") as fh:
        b = fh.read()
    dst = os.path.join(HW, "calibration", f"ibm_kingston_properties_{CALIBRATION_STAMP}.json.gz")
    out = write_gz_json(dst, to_utc(json.loads(b)))
    manifest.append({"path": rel(dst), "what": "device calibration record in force when the 4 jobs ran",
                     "source": "Stage8_Hardware/shared/snapshots/ibm_kingston/jobs/props_<job id>.json (4 identical files)",
                     "source_sha256": hashes.pop(), "sha256": sha256_bytes(out), "bytes": len(out),
                     "transform": "time stamps converted to UTC; gzip"})

    # 3. frozen circuits, byte for byte
    for jid, (role, name) in JOBS.items():
        for ext in (".json", ".qpy", "_params.npz"):
            s_rel = f"{FROZEN_DIR[name[:2]]}/{name}{ext}"
            dst = os.path.join(HW, "circuits", f"{name}{ext}")
            os.makedirs(os.path.dirname(dst), exist_ok=True)
            shutil.copyfile(src(s_rel), dst)
            manifest.append({"path": rel(dst), "what": f"frozen job {name} ({role}): " +
                             {".json": "label index, shots, hash", ".qpy": "transpiled circuits (QPY)",
                              "_params.npz": "parameter values of every binding"}[ext],
                             "source": s_rel, "source_sha256": sha256_file(src(s_rel)), "sha256": sha256_file(dst),
                             "bytes": os.path.getsize(dst), "transform": "none (byte copy)"})

    # 4. pre-registration record, byte for byte; ledger reduced to the 4 hardware jobs
    for s_rel, name in PREREG:
        dst = os.path.join(HW, "preregistration", name)
        os.makedirs(os.path.dirname(dst), exist_ok=True)
        shutil.copyfile(src(s_rel), dst)
        manifest.append({"path": rel(dst), "what": "pre-registration record", "source": s_rel,
                         "source_sha256": sha256_file(src(s_rel)), "sha256": sha256_file(dst),
                         "bytes": os.path.getsize(dst), "transform": "none (byte copy)"})
    s_rel = "Stage8_Hardware/shared/jobs/hw_jobs.json"
    with open(src(s_rel), "rb") as fh:
        b = fh.read()
    jobs = [{k: j[k] for k in RECORD_KEEP + ("status",) if k in j}
            for j in json.loads(b)["jobs"] if j["where"] == "hardware"]
    assert [j["id"] for j in jobs] == list(JOBS)
    dst = os.path.join(HW, "preregistration", "ledger.json")
    with open(dst, "w", encoding="utf8") as fh:
        json.dump({"jobs": jobs}, fh, indent=1)
    manifest.append({"path": rel(dst), "what": "job ledger: the 4 hardware jobs", "source": s_rel,
                     "source_sha256": sha256_bytes(b), "sha256": sha256_file(dst), "bytes": os.path.getsize(dst),
                     "transform": "local simulator runs removed; experiment fields only"})

    # 5. reference results
    for s_rel, name, cleaned in REFERENCE:
        dst = os.path.join(ref_dir, name)
        os.makedirs(ref_dir, exist_ok=True)
        if cleaned:
            with open(src(s_rel), encoding="utf8") as fh:
                obj = json.load(fh)
            if cleaned == "numbers":
                obj, how = clean_numbers(obj), "the plan's quota entries (limit, remaining) removed"
            else:
                obj, how = clean_reference(obj), "account quota fields and the usage-detail metric id removed"
            with open(dst, "w", encoding="utf8") as fh:
                json.dump(obj, fh, indent=1, ensure_ascii=False)
        else:
            shutil.copyfile(src(s_rel), dst)
            how = "none (byte copy)"
        manifest.append({"path": rel(dst), "what": "saved result used as a regression target", "source": s_rel,
                         "source_sha256": sha256_file(src(s_rel)), "sha256": sha256_file(dst),
                         "bytes": os.path.getsize(dst), "transform": how})

    # 6. saved simulation results, byte for byte
    for s_rel, name in SIMULATION:
        dst = os.path.join(ref_dir, "simulation", name)
        os.makedirs(os.path.dirname(dst), exist_ok=True)
        shutil.copyfile(src(s_rel), dst)
        manifest.append({"path": rel(dst), "what": "saved simulation result used by the notebooks", "source": s_rel,
                         "source_sha256": sha256_file(src(s_rel)), "sha256": sha256_file(dst),
                         "bytes": os.path.getsize(dst), "transform": "none (byte copy)"})

    # 7. check: cleaned counts == original counts, pub by pub
    for jid in JOBS:
        with gzip.open(os.path.join(HW, "counts", f"{jid}.json.gz"), "rt", encoding="utf8") as fh:
            c = json.load(fh)
        with open(src(f"Stage8_Hardware/shared/jobs/raw/{jid}.json"), encoding="utf8") as fh:
            o = json.load(fh)
        assert c["pubs"] == o["pubs"]
    with open(os.path.join(DATA, "MANIFEST.json"), "w", encoding="utf8") as fh:
        json.dump({"note": "source paths are relative to the research project folder (not shipped)",
                   "files": manifest}, fh, indent=1)
    print(f"{len(manifest)} files written; counts identical to the originals for all {len(JOBS)} jobs")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--project", required=True, help="path of the research project folder")
    main(os.path.abspath(ap.parse_args().project))
