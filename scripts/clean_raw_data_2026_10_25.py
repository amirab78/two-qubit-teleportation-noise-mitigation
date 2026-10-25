"""Build ``data/hardware/ibm_kingston_ibm_fez_2026-10-25/`` and ``data/reference/replication/`` from the research
project's extra hardware stage (run once by the authors; kept to document how the shipped data were made).

usage:  python scripts/clean_raw_data_2026_10_25.py --project PATH_TO_RESEARCH_PROJECT

The same rules as ``clean_raw_data.py`` (the 4 October data), applied to the 11 jobs of 25 October 2026 on
``ibm_kingston`` and ``ibm_fez``:

* counts -> ``counts/<job id>.json.gz``: counts and per-shot outcomes unchanged; record and metrics reduced to the
  fields that describe the experiment (the usage-detail metric id is dropped);
* the job-time calibration records (one per distinct record: 4 on Kingston, 1 on Fez) -> gzip, time stamps
  converted from the local time zone of the recording computer to UTC;
* frozen circuits, predictions, their SHA-256 record and the approvals: byte copies, so every recorded hash can be
  checked; the job ledger reduced to experiment fields;
* the saved per-job analyses ``x_<job>.json`` -> ``data/reference/replication/`` with the IBM plan's quota fields
  (``budget``, ``budget_now``) removed;
* ``data/MANIFEST.json``: the entries of this run are replaced (the 4 October entries are left as they are).

Only the folders of this run are written; the 4 October data are not touched.
"""
from __future__ import annotations

import argparse
import gzip
import json
import os
import shutil

from clean_raw_data import DATA, REPO, sha256_bytes, sha256_file, to_utc, write_gz_json

RUN = "ibm_kingston_ibm_fez_2026-10-25"
HW = os.path.join(DATA, "hardware", RUN)
REF = os.path.join(DATA, "reference", "replication")
SRC = "Extra_Simulations_Real_Quantum_Computer"

JOBS = {  # job id -> frozen job name
    "dbeqnffsstms73fciu80": "K0_ibm_kingston", "dbeqouvgdj9s73fma3c0": "K1_ibm_kingston",
    "dber02gm94ac738sb5tg": "K2_ibm_kingston", "dber2sfpc3tc739m7j00": "K3_ibm_kingston",
    "dbere1nsstms73fck4pg": "K4_ibm_kingston", "dberne7gdj9s73fmbpd0": "F0_ibm_fez",
    "dberor7pc3tc739m8tlg": "F1_ibm_fez", "dberrs8m94ac738scpig": "F2_ibm_fez",
    "dbes408m94ac738sd89g": "G1_ibm_fez", "dbes5lngdj9s73fmcl7g": "G2_ibm_kingston",
    "dbes9k7pc3tc739m9vkg": "G3_ibm_kingston",
}
RECORD_KEEP = ("id", "where", "backend", "round", "purpose", "name", "labels", "n_pubs", "shots", "shots_total",
               "estimate_s", "max_seconds", "dd", "twirl", "x_options", "sampler", "submitted_utc", "usage_s", "hash",
               "xhash", "ibm_estimate_s")
QUOTA_KEYS = {"budget", "budget_now"}       # account quota of the IBM plan: not part of the experiment
PREDICTIONS = ("x_pred_K0K1.json", "x_pred_K234_ibm_kingston.json", "x_pred_F0F1_ibm_fez.json",
               "x_pred_F2_ibm_fez.json", "x_pred_G.json", "x_pred_G3.json", "sha256.json")
RESULTS = ["x_K0_ibm_kingston.json", "x_K1_ibm_kingston.json", "x_K2_ibm_kingston.json", "x_K3_ibm_kingston.json",
           "x_K4_ibm_kingston.json", "x_F0_ibm_fez.json", "x_F1_ibm_fez.json", "x_F2_ibm_fez.json",
           "x_G1_ibm_fez.json", "x_G2_ibm_kingston.json", "x_G3_ibm_kingston.json", "x_K1_retro_wait.json",
           "x_l1_tests_ibm_kingston.json", "x_l1_tests_ibm_fez.json", "x_l1_tests_G.json", "x_l1_tests_G3.json"]


def clean_metrics(m: dict | None) -> dict | None:
    if not m:
        return m
    out = {k: m[k] for k in ("timestamps", "circuits_execution_time_ns", "qiskit_version") if k in m}
    if "usage" in m:
        out["usage"] = {k: v for k, v in m["usage"].items() if k != "details"}
    return out


def drop_quota(obj):
    if isinstance(obj, dict):
        return {k: drop_quota(v) for k, v in obj.items() if k not in QUOTA_KEYS}
    if isinstance(obj, list):
        return [drop_quota(v) for v in obj]
    return obj


def main(project: str) -> None:
    src = lambda rel: os.path.join(project, SRC, *rel.split("/"))        # noqa: E731
    rel = lambda p: os.path.relpath(p, REPO).replace(os.sep, "/")        # noqa: E731
    srel = lambda r: f"{SRC}/{r}"                                        # noqa: E731
    manifest = []
    for d in (HW, REF):
        if os.path.isdir(d):
            shutil.rmtree(d)

    def add(dst, what, source, source_sha, how):
        manifest.append({"path": rel(dst), "what": what, "source": source, "source_sha256": source_sha,
                         "sha256": sha256_file(dst), "bytes": os.path.getsize(dst), "transform": how})

    def copy(s_rel, dst, what):
        os.makedirs(os.path.dirname(dst), exist_ok=True)
        shutil.copyfile(src(s_rel), dst)
        add(dst, what, srel(s_rel), sha256_file(src(s_rel)), "none (byte copy)")

    # 1. counts of the 11 hardware jobs
    with open(src("jobs/hw_jobs.json"), "rb") as fh:
        ledger_bytes = fh.read()
    ledger = {j["id"]: j for j in json.loads(ledger_bytes)["jobs"] if j["where"] == "hardware"}
    assert list(ledger) == list(JOBS), "the ledger must list exactly the 11 jobs, in order"
    for jid, name in JOBS.items():
        with open(src(f"jobs/raw/{jid}.json"), "rb") as fh:
            b = fh.read()
        raw = json.loads(b)
        rec = dict(raw["record"], usage_s=ledger[jid]["usage_s"])        # the raw record was saved before usage was known
        cleaned = {"id": raw["id"], "record": {k: rec[k] for k in RECORD_KEEP if k in rec}, "pubs": raw["pubs"],
                   "metrics": clean_metrics(raw.get("metrics"))}
        assert cleaned["pubs"] == raw["pubs"], "counts must be unchanged"
        dst = os.path.join(HW, "counts", f"{jid}.json.gz")
        write_gz_json(dst, cleaned)
        add(dst, f"counts and per-shot outcomes, {name}", srel(f"jobs/raw/{jid}.json"), sha256_bytes(b),
            "kept counts unchanged; record and metrics reduced to experiment fields; usage from the ledger; gzip")

    # 2. the job-time calibration records (one file per distinct record)
    seen: dict[str, str] = {}
    for jid, name in JOBS.items():
        dev = name.split("_", 1)[1]
        p = src(f"snapshots/{dev}/jobs/props_{jid}.json")
        h = sha256_file(p)
        if h in seen:
            continue
        stamp = next(s for s in sorted(os.listdir(src(f"snapshots/{dev}")))
                     if s != "jobs" and sha256_file(src(f"snapshots/{dev}/{s}/props.json")) == h)
        seen[h] = stamp
        with open(p, "rb") as fh:
            b = fh.read()
        dst = os.path.join(HW, "calibration", f"{dev}_properties_{stamp}.json.gz")
        write_gz_json(dst, to_utc(json.loads(b)))
        users = [n for j, n in JOBS.items() if sha256_file(src(f"snapshots/{n.split('_', 1)[1]}/jobs/props_{j}.json")) == h]
        add(dst, f"device calibration record in force for {', '.join(users)}",
            srel(f"snapshots/{dev}/jobs/props_<job id>.json"), h, "time stamps converted to UTC; gzip")

    # 3. frozen circuits, byte for byte
    for name in JOBS.values():
        for ext, what in ((".json", "label index, shots, hash, prediction file and its SHA-256"),
                          (".qpy", "transpiled circuits (QPY)"), ("_params.npz", "parameter values of every binding")):
            copy(f"frozen/{name}{ext}", os.path.join(HW, "circuits", f"{name}{ext}"), f"frozen job {name}: {what}")

    # 4. pre-registration: predictions, their SHA-256 record, approvals (byte copies); the ledger (experiment fields)
    for name in PREDICTIONS:
        copy(f"predictions/{name}", os.path.join(HW, "preregistration", name), "pre-registration record")
    copy("jobs/approved.json", os.path.join(HW, "preregistration", "approvals.json"), "approvals (hash and time cap)")
    jobs = [{k: j[k] for k in RECORD_KEEP + ("status",) if k in j} for j in ledger.values()]
    dst = os.path.join(HW, "preregistration", "ledger.json")
    with open(dst, "w", encoding="utf8") as fh:
        json.dump({"jobs": jobs}, fh, indent=1)
    add(dst, "job ledger: the 11 hardware jobs", srel("jobs/hw_jobs.json"), sha256_bytes(ledger_bytes),
        "experiment fields only")

    # 5. the saved per-job analyses (regression targets), quota fields removed
    os.makedirs(REF, exist_ok=True)
    for name in RESULTS:
        with open(src(f"results/{name}"), encoding="utf8") as fh:
            obj = json.load(fh)
        dst = os.path.join(REF, name)
        with open(dst, "w", encoding="utf8") as fh:
            json.dump(drop_quota(obj), fh, indent=1, ensure_ascii=False)
        add(dst, "saved result used as a regression target", srel(f"results/{name}"),
            sha256_file(src(f"results/{name}")), "the IBM plan's quota fields (budget, budget_now) removed")

    # 6. check: cleaned counts == original counts, pub by pub
    for jid in JOBS:
        with gzip.open(os.path.join(HW, "counts", f"{jid}.json.gz"), "rt", encoding="utf8") as fh:
            c = json.load(fh)
        with open(src(f"jobs/raw/{jid}.json"), encoding="utf8") as fh:
            assert c["pubs"] == json.load(fh)["pubs"]

    path = os.path.join(DATA, "MANIFEST.json")
    with open(path, encoding="utf8") as fh:
        man = json.load(fh)
    keep = [e for e in man["files"] if not e["path"].startswith((rel(HW) + "/", rel(REF) + "/"))]
    man["files"] = keep + manifest
    with open(path, "w", encoding="utf8") as fh:
        json.dump(man, fh, indent=1)
    print(f"{len(manifest)} files written; counts identical to the originals for all {len(JOBS)} jobs")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--project", required=True, help="path of the research project folder")
    main(os.path.abspath(ap.parse_args().project))
