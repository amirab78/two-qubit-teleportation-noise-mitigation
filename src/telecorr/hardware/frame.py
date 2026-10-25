"""The Pauli frame: Bob's corrections applied in software (circuit version A).

Roadmap: step 8 (real hardware).

Alice's register (4 bits, bit 0 = q0, bit 1 = q1, bit 2 = q2, bit 3 = q3) gives the Z corrections of Bob's q4
(bit 0) and q5 (bit 1) and the X corrections of q4 (bit 2) and q5 (bit 3). Measuring Bob in basis B, an X
correction flips Z and Y results, a Z correction flips X and Y results:

    flip(Z) = x,  flip(X) = z,  flip(Y) = x XOR z.
"""
from __future__ import annotations

import numpy as np


def frame_bit(bob_bit, basis: str, z_bit, x_bit):
    """Corrected value of one of Bob's bits measured in ``basis``."""
    flip = {"Z": x_bit, "X": z_bit, "Y": x_bit ^ z_bit}[basis]
    return bob_bit ^ flip


def frame_keys(alice: np.ndarray, bob: np.ndarray, bases: str) -> np.ndarray:
    """Corrected Bob outcome (integer b5 * 2 + b4) of every shot; ``bases`` = (q4 basis, q5 basis)."""
    a = np.asarray(alice); b = np.asarray(bob)
    z4, z5, x4, x5 = a & 1, (a >> 1) & 1, (a >> 2) & 1, (a >> 3) & 1
    b4 = frame_bit(b & 1, bases[0], z4, x4)
    b5 = frame_bit((b >> 1) & 1, bases[1], z5, x5)
    return b5 * 2 + b4


def frame_counts(alice: np.ndarray, bob: np.ndarray, bases: str) -> dict[str, int]:
    """Counts {"b5b4": n} of Bob's frame-corrected outcomes."""
    c = np.bincount(frame_keys(alice, bob, bases), minlength=4)
    return {f"{k >> 1}{k & 1}": int(c[k]) for k in range(4)}


def raw_counts(bob: np.ndarray) -> dict[str, int]:
    """Counts of Bob's bits as measured (versions B and C correct physically)."""
    c = np.bincount(np.asarray(bob), minlength=4)
    return {f"{k >> 1}{k & 1}": int(c[k]) for k in range(4)}
