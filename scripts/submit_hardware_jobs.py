"""How the hardware jobs of the study were sent. DOCUMENTATION ONLY: disabled by default.

Nothing in this repository needs a quantum computer; every hardware result is rebuilt from the shipped counts
(``python -m telecorr.hardware reproduce``). This script shows the send procedure of the study, so that it can be
reviewed or repeated on another device with a new, separately approved job:

1. load a frozen job (``data/hardware/<run>/circuits/<name>.qpy`` + ``_params.npz`` + ``.json`` manifest);
2. check the manifest against the approvals: the job's hash must be approved, with exactly the manifest's time cap,
   and it must never have been sent before (the job ledger);
3. set every sampler option explicitly (no server defaults): the time cap (``max_execution_time``), dynamical
   decoupling and twirling as frozen;
4. send the PUBs (circuit, parameter values, shots) with SamplerV2 and record the job id.

The study's guard additionally recomputed the job hash from the reloaded circuits (a canonical text of every
instruction); that step is not repeated here.

The script refuses to run unless BOTH ``--enable-hardware`` is given AND the environment variable
``TELECORR_ENABLE_HARDWARE=1`` is set. Credentials are never read from this file or from the command line: the
IBM Quantum account must already be saved on the local machine (``QiskitRuntimeService.save_account``), and only its
local *name* is passed with ``--account``. Requires the optional package ``qiskit-ibm-runtime`` (not installed by
this repository).

usage (documentation):  python scripts/submit_hardware_jobs.py --job R1_main_ibm_kingston --account <saved name>
                         --enable-hardware          (and TELECORR_ENABLE_HARDWARE=1 in the environment)
"""
from __future__ import annotations

import argparse
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
RUN = os.path.join(os.path.dirname(HERE), "data", "hardware", "ibm_kingston_2026-10-04")


def load_frozen(name: str):
    """The frozen job: PUBs (circuit, parameter values or None, shots) and its manifest."""
    import numpy as np
    from qiskit import qpy
    base = os.path.join(RUN, "circuits", name)
    with open(base + ".json", encoding="utf8") as fh:
        man = json.load(fh)
    with open(base + ".qpy", "rb") as fh:
        circuits = qpy.load(fh)
    arrs = np.load(base + "_params.npz")
    pubs = [(c, arrs[f"p{i}"] if has else None, shots)
            for i, (c, has, shots) in enumerate(zip(circuits, man["has_params"], man["shots"]))]
    if len(pubs) != man["n_pubs"]:
        raise ValueError(f"{name}: {len(pubs)} PUBs reloaded, the manifest lists {man['n_pubs']}")
    return pubs, man


def check_approved(man: dict) -> None:
    """The send guard: approved by hash, with exactly this time cap, and never sent before."""
    with open(os.path.join(RUN, "preregistration", "approvals.json"), encoding="utf8") as fh:
        approved = [a for a in json.load(fh)["approved"] if a["hash"] == man["hash"]]
    if not approved:
        raise PermissionError(f"job {man['name']} (hash {man['hash'][:12]}) is not approved; nothing was sent")
    if int(man["max_seconds"]) != int(approved[-1]["max_seconds"]):
        raise PermissionError(f"job {man['name']}: time cap {man['max_seconds']} s != approved {approved[-1]['max_seconds']} s")
    with open(os.path.join(RUN, "preregistration", "ledger.json"), encoding="utf8") as fh:
        if any(j.get("hash") == man["hash"] for j in json.load(fh)["jobs"]):
            raise PermissionError(f"job {man['name']} was already sent once (see ledger.json); nothing was sent")


def send(pubs, man: dict, account: str) -> str:
    from qiskit_ibm_runtime import QiskitRuntimeService, SamplerV2          # optional dependency
    service = QiskitRuntimeService(name=account)                            # a locally saved account, by name only
    backend = service.backend(man["device"])
    sampler = SamplerV2(mode=backend)
    sampler.options.max_execution_time = int(man["max_seconds"])
    sampler.options.dynamical_decoupling.enable = bool(man["options"]["dd"])
    sampler.options.twirling.enable_gates = bool(man["options"]["twirl"])
    sampler.options.twirling.enable_measure = bool(man["options"]["twirl"])
    job = sampler.run(pubs)
    return job.job_id()


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Send a frozen, approved hardware job (documentation only; disabled by default).")
    ap.add_argument("--job", required=True, help="frozen job name, e.g. R1_main_ibm_kingston")
    ap.add_argument("--account", required=True, help="NAME of an IBM Quantum account saved on this machine (no key)")
    ap.add_argument("--enable-hardware", action="store_true", help="required, together with TELECORR_ENABLE_HARDWARE=1")
    a = ap.parse_args(argv)
    if not (a.enable_hardware and os.environ.get("TELECORR_ENABLE_HARDWARE") == "1"):
        print("Disabled: this script documents how the study's jobs were sent and does not run by default.\n"
              "Every hardware result is rebuilt from the shipped counts: python -m telecorr.hardware reproduce")
        return 2
    pubs, man = load_frozen(a.job)
    check_approved(man)
    print(f"sending {man['name']} to {man['device']}: {man['n_pubs']} PUBs, {man['shots_total']} shots, cap {man['max_seconds']} s")
    print("job id:", send(pubs, man, a.account))
    return 0


if __name__ == "__main__":
    sys.exit(main())
