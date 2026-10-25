"""The shipped data are what the research project recorded: file hashes, the pre-registration chain
(predictions hashed before the jobs were sent; the job hash approved = frozen = sent), and no account information."""
import gzip
import hashlib
import json
import re

from telecorr.hardware import data as D


def sha(path):
    with open(path, "rb") as fh:
        return hashlib.sha256(fh.read()).hexdigest()


def test_manifest_hashes():
    root = D.data_dir().parent
    man = json.load(open(D.data_dir() / "MANIFEST.json", encoding="utf8"))
    assert len(man["files"]) >= 40
    for f in man["files"]:
        assert sha(root / f["path"]) == f["sha256"], f["path"]
        if f["transform"].startswith("none"):
            assert f["sha256"] == f["source_sha256"], f["path"]       # byte copies keep the source hash


def test_predictions_hashed_before_sending():
    pre = D.run_dir() / "preregistration"
    for record in ("8c_predictions_sha256.json", "8d_predictions_sha256.json"):
        rec = json.load(open(pre / record, encoding="utf8"))
        names = [k for k in rec if k.endswith(".json")]
        assert names
        for name in names:
            assert sha(pre / name) == rec[name], name
    t8c = json.load(open(pre / "8c_predictions_sha256.json"))["written_utc"]
    t8d = json.load(open(pre / "8d_predictions_sha256.json"))
    ledger = {j["purpose"]: j for j in json.load(open(pre / "ledger.json"))["jobs"]}
    assert t8c < ledger["8d_round1"]["submitted_utc"]                 # 8c predictions before Round 1
    assert t8d["written_utc"] < ledger["8d_round2_main"]["submitted_utc"]
    assert t8d["tomo_model_written_utc"] < ledger["8d_round2_main"]["submitted_utc"]


def test_job_hash_chain():
    """approved hash = frozen job hash = hash recorded when the job was sent, for all 4 jobs."""
    approved = {a["name"]: a for a in D.preregistration("approvals.json")["approved"]}
    for purpose, (jid, name) in D.JOBS.items():
        frozen = D.frozen_job(purpose)
        raw = D.raw_job(purpose)
        assert raw["id"] == jid
        assert approved[name]["hash"] == frozen["hash"] == raw["record"]["hash"], name
        assert raw["record"]["max_seconds"] == approved[name]["max_seconds"]
        assert raw["record"]["usage_s"] <= raw["record"]["max_seconds"]


def test_counts_are_complete():
    for purpose in D.JOBS:
        raw = D.raw_job(purpose); idx = D.frozen_job(purpose)["index"]
        assert [p["label"] for p in raw["pubs"]] == list(idx)
        for p in raw["pubs"]:
            shots = int(p["label"].split("@")[1])
            for reg, c in p["counts"].items():
                tables = c if isinstance(c, list) else [c]
                assert len(tables) == len(idx[p["label"]])
                assert all(sum(t.values()) == shots for t in tables)


def test_no_account_information_in_data():
    pat = re.compile(r"(crn:|apikey|api_key|token|instance|@[a-z0-9-]+\.[a-z]|usage_limit|usage_period|"
                     r"[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}|\+03:30)", re.I)
    for path in D.data_dir().rglob("*"):
        if path.suffix == ".gz":
            text = gzip.open(path, "rt", encoding="utf8").read()
        elif path.suffix == ".json":
            text = path.read_text(encoding="utf8")
        else:
            continue
        hits = {m.group(0) for m in pat.finditer(text)}
        assert not hits, (path.name, sorted(hits)[:5])


# ------------------------------------------------------------------------------- the replication of 25 October 2026
def test_replication_predictions_hashed_before_sending():
    """Every prediction file of the replication has the SHA-256 recorded before sending, the frozen job names that
    file and hash, and the file was written before the job was submitted."""
    pre = D.replication_dir() / "preregistration"
    rec = json.load(open(pre / "sha256.json", encoding="utf8"))
    for name, r in rec.items():
        assert sha(pre / name) == r["sha256"], name
    for job in D.REPLICATION_JOBS:
        man = D.replication_frozen(job); raw = D.replication_raw(job)
        assert sha(pre / man["prediction_file"]) == man["prediction_sha256"], job
        assert rec[man["prediction_file"]]["written_utc"] < raw["record"]["submitted_utc"], job


def test_replication_job_hash_chain():
    """approved hash = frozen job hash = hash recorded when the job was sent, and usage <= the approved cap."""
    approved = {a["name"]: a for a in D.replication_preregistration("approvals.json")["approved"]}
    ledger = {j["id"]: j for j in D.replication_preregistration("ledger.json")["jobs"]}
    for job, (jid, dev) in D.REPLICATION_JOBS.items():
        name = f"{job}_{dev}"
        man = D.replication_frozen(job); raw = D.replication_raw(job)
        assert raw["id"] == jid and ledger[jid]["name"] == name
        assert approved[name]["xhash"] == man["xhash"] == raw["record"]["xhash"] == ledger[jid]["xhash"], name
        assert man["hash"] == raw["record"]["hash"], name
        assert raw["record"]["max_seconds"] == approved[name]["max_seconds"] == man["max_seconds"]
        assert raw["record"]["usage_s"] <= raw["record"]["max_seconds"]
        assert approved[name]["approved_utc"] < raw["record"]["submitted_utc"]


def test_replication_counts_are_complete():
    for job in D.REPLICATION_JOBS:
        raw = D.replication_raw(job); idx = D.replication_frozen(job)["index"]
        assert [p["label"] for p in raw["pubs"]] == list(idx)
        for p in raw["pubs"]:
            shots = int(p["label"].split("@")[1])
            for reg, c in p["counts"].items():
                tables = c if isinstance(c, list) else [c]
                assert len(tables) == len(idx[p["label"]])
                assert all(sum(t.values()) == shots for t in tables)
