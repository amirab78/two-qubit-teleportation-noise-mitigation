"""Loading the shipped hardware data (``data/hardware/ibm_kingston_2026-10-04``).

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
