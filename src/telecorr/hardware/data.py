"""Loading the shipped hardware data (``data/hardware/ibm_kingston_2026-10-04`` and the replication of
25 October 2026, ``data/hardware/ibm_kingston_ibm_fez_2026-10-25``).

Roadmap: step 8 (real hardware).

The data folder is found from the environment variable ``TELECORR_DATA`` (path of the repository's ``data``
folder) or, for a source checkout, next to ``src/``.
"""
from __future__ import annotations

import gzip
import json
import os
from functools import lru_cache
from pathlib import Path

DEVICE = "ibm_kingston"
RUN = "ibm_kingston_2026-10-04"
#: purpose -> (job id, frozen job name)
JOBS: dict[str, tuple[str, str]] = {
    "pilot": ("db16kn2vog1s73fhnpog", "R0_main_ibm_kingston"),
    "round1": ("db16p1bid5ic73eqnsug", "R1_main_ibm_kingston"),
    "round2_main": ("db175bjid5ic73eqoav0", "R2_main_ibm_kingston"),
    "round2_dd": ("db175givog1s73fhocjg", "R2_dd_ibm_kingston"),
}


def data_dir() -> Path:
    env = os.environ.get("TELECORR_DATA")
    p = Path(env) if env else Path(__file__).resolve().parents[3] / "data"
    if not (p / "hardware" / RUN).is_dir():
        raise FileNotFoundError(f"hardware data not found under {p}; set TELECORR_DATA to the repository's data folder")
    return p


def run_dir() -> Path:
    return data_dir() / "hardware" / RUN


def reference_dir() -> Path:
    return data_dir() / "reference"


def _json(path: Path):
    if path.suffix == ".gz":
        with gzip.open(path, "rt", encoding="utf8") as fh:
            return json.load(fh)
    with open(path, encoding="utf8") as fh:
        return json.load(fh)


@lru_cache(maxsize=None)
def raw_job(purpose: str) -> dict:
    """The cleaned raw result of one hardware job: {"id", "record", "pubs", "metrics"}."""
    jid, _ = JOBS[purpose]
    return _json(run_dir() / "counts" / f"{jid}.json.gz")


@lru_cache(maxsize=None)
def frozen_job(purpose: str) -> dict:
    """The frozen job manifest (labels, label index, shots, hash) written before the job was sent."""
    _, name = JOBS[purpose]
    return _json(run_dir() / "circuits" / f"{name}.json")


def preregistration(name: str) -> dict:
    """A pre-registration file, e.g. ``predictions_round1.json`` or ``approvals.json``."""
    return _json(run_dir() / "preregistration" / name)


def reference(name: str) -> dict:
    """A saved result of the research project used as a regression target (``data/reference``)."""
    return _json(reference_dir() / name)


def calibration() -> dict:
    """The device calibration record in force when the jobs ran (time stamps in UTC)."""
    files = sorted((run_dir() / "calibration").glob("*.json.gz"))
    return _json(files[0])


def round_data(raw: dict, index: dict) -> dict:
    """{binding label: {"ints": {register: per-shot integers}, "counts": {register: counts}}} of one job.
    ``index`` maps each pub label to the binding labels of its parameter sets (from the frozen manifest)."""
    out = {}
    for p in raw["pubs"]:
        labs = index[p["label"]]
        regs = list(p["counts"])
        for i, lab in enumerate(labs):
            if isinstance(p["counts"][regs[0]], list):
                out[lab] = {"ints": {r: p["shots_int"][r][i] for r in regs}, "counts": {r: p["counts"][r][i] for r in regs}}
            else:
                out[lab] = {"ints": {r: p["shots_int"][r] for r in regs}, "counts": {r: p["counts"][r] for r in regs}}
    return out


def pub_of_label(purpose: str) -> dict[str, str]:
    """binding label -> pub label (``"A_product@3000"`` etc.) for one job."""
    idx = frozen_job(purpose)["index"]
    return {lab: pub for pub, labs in idx.items() for lab in labs}


def load_round(n: int) -> dict:
    """Binding-level data of Round 1 (``n = 1``) or Round 2 (``n = 2``: main job, then the DD job)."""
    purposes = {1: ("round1",), 2: ("round2_main", "round2_dd")}[n]
    out: dict = {}
    for p in purposes:
        out.update(round_data(raw_job(p), frozen_job(p)["index"]))
    return out


def _version_of(pub: str) -> str:
    isa = pub.split("@")[0]
    return "msg" if isa.startswith("alice_") else ("A" if isa.startswith("A_") else "raw")


# --------------------------------------------------------------------------------------- the replication (25 Oct)
REPLICATION_RUN = "ibm_kingston_ibm_fez_2026-10-25"
#: job -> IBM job id; the frozen job is "<job>_<device>"
REPLICATION_JOBS: dict[str, tuple[str, str]] = {
    "K0": ("dbeqnffsstms73fciu80", "ibm_kingston"), "K1": ("dbeqouvgdj9s73fma3c0", "ibm_kingston"),
    "K2": ("dber02gm94ac738sb5tg", "ibm_kingston"), "K3": ("dber2sfpc3tc739m7j00", "ibm_kingston"),
    "K4": ("dbere1nsstms73fck4pg", "ibm_kingston"), "F0": ("dberne7gdj9s73fmbpd0", "ibm_fez"),
    "F1": ("dberor7pc3tc739m8tlg", "ibm_fez"), "F2": ("dberrs8m94ac738scpig", "ibm_fez"),
    "G1": ("dbes408m94ac738sd89g", "ibm_fez"), "G2": ("dbes5lngdj9s73fmcl7g", "ibm_kingston"),
    "G3": ("dbes9k7pc3tc739m9vkg", "ibm_kingston"),
}


def replication_dir() -> Path:
    return data_dir() / "hardware" / REPLICATION_RUN


@lru_cache(maxsize=None)
def replication_raw(job: str) -> dict:
    """The cleaned raw result of one replication job (``"K1"``, ``"F2"``, ``"G3"`` ...)."""
    return _json(replication_dir() / "counts" / f"{REPLICATION_JOBS[job][0]}.json.gz")


@lru_cache(maxsize=None)
def replication_frozen(job: str) -> dict:
    """The frozen job manifest (labels, label index, shots, hash, prediction file and its SHA-256)."""
    return _json(replication_dir() / "circuits" / f"{job}_{REPLICATION_JOBS[job][1]}.json")


def replication_preregistration(name: str) -> dict:
    """A pre-registration file of the replication, e.g. ``x_pred_K0K1.json`` or ``approvals.json``."""
    return _json(replication_dir() / "preregistration" / name)


def replication_reference(name: str) -> dict:
    """A saved per-job analysis of the replication (``data/reference/replication/x_<job>.json``)."""
    return _json(reference_dir() / "replication" / name)


def replication_calibration(device: str, stamp: str) -> dict:
    """A job-time calibration record of the replication (time stamps in UTC)."""
    return _json(replication_dir() / "calibration" / f"{device}_properties_{stamp}.json.gz")


def replication_job(job: str) -> tuple[dict, dict[str, int], dict[str, str]]:
    """Binding-level data of one replication job, the shots of every binding and its post-processing version."""
    man = replication_frozen(job)
    data = round_data(replication_raw(job), man["index"])
    shots, versions = {}, {}
    for pub, labs in man["index"].items():
        for lab in labs:
            shots[lab] = int(pub.split("@")[1]); versions[lab] = _version_of(pub)
    return data, shots, versions


def labels_info(n: int) -> tuple[dict[str, str], dict[str, int]]:
    """For every binding label of a round: its post-processing version ("A" = Pauli frame, "raw" = Bob's bits
    as measured in versions B / C and the readout circuits, "msg" = Alice reference) and its shots."""
    purposes = {1: ("round1",), 2: ("round2_main", "round2_dd")}[n]
    versions, shots = {}, {}
    for p in purposes:
        for lab, pub in pub_of_label(p).items():
            isa, n_shots = pub.split("@")
            versions[lab] = "msg" if isa.startswith("alice_") else ("A" if isa.startswith("A_") else "raw")
            shots[lab] = int(n_shots)
    return versions, shots
