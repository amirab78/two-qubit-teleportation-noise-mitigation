"""Readout calibration from the four pattern circuits and Bob-only readout mitigation.

Roadmap: step 8 (real hardware).
"""
from __future__ import annotations

from typing import Mapping

import numpy as np

from .frame import frame_keys

READOUT_PATTERNS: tuple[str, ...] = ("000000", "111111", "010101", "101010")   # logical q0..q5, left = q0


def confusion_from_counts(counts_by_pattern: Mapping[str, Mapping[str, Mapping[str, int]]]) -> dict:
    """Per-qubit e0 = P(1|0) and e1 = P(0|1) (logical q0..q5) and Bob's 4 x 4 matrix M[measured, prepared]
    (index b5 * 2 + b4). ``counts_by_pattern[pattern] = {"alice": counts, "bob": counts}``."""
    e = {i: {"n0": 0, "f0": 0, "n1": 0, "f1": 0} for i in range(6)}
    M = np.zeros((4, 4))
    for pat, cc in counts_by_pattern.items():
        for reg, off, nb in (("alice", 0, 4), ("bob", 4, 2)):
            for key, cnt in cc[reg].items():
                for j in range(nb):
                    i = off + j; prep = int(pat[i]); got = int(key[-1 - j])
                    e[i][f"n{prep}"] += cnt; e[i][f"f{prep}"] += cnt * (got != prep)
        col = int(pat[5]) * 2 + int(pat[4])
        tot = sum(cc["bob"].values())
        for key, cnt in cc["bob"].items():
            M[int(key, 2), col] += cnt / tot
    return {"e0": [e[i]["f0"] / e[i]["n0"] for i in range(6)], "e1": [e[i]["f1"] / e[i]["n1"] for i in range(6)],
            "M_bob": M.tolist()}


def mitigate_bob(alice: np.ndarray, bob: np.ndarray, bases: str, M, version: str = "A") -> dict[str, float]:
    """Bob-only readout mitigation, returning a quasi-probability {"b5b4": p} (entries may be < 0).

    The confusion matrix acts on Bob's RAW bits, so in version A it is inverted for each of Alice's outcomes
    separately and only then is the Pauli frame applied (exact). Alice's readout errors are a real teleportation
    error and are not mitigated."""
    Minv = np.linalg.inv(np.asarray(M))
    a = np.asarray(alice); b = np.asarray(bob)
    out = np.zeros(4)
    for m in (np.unique(a) if version == "A" else [None]):
        sel = (a == m) if m is not None else np.ones_like(b, bool)
        q = Minv @ np.bincount(b[sel], minlength=4).astype(float)
        if version == "A":
            np.add.at(out, frame_keys(np.full(4, m), np.arange(4), bases), q)
        else:
            out += q
    out /= out.sum()
    return {f"{k >> 1}{k & 1}": float(out[k]) for k in range(4)}
