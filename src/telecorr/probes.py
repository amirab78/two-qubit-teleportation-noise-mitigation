"""The counts-only probe cost: what the calibrator measures instead of the (unknown) fidelity.

Roadmap: step 3 (correction from counts); the same circuits are TOMO-9's in steps 7 and 8.

Three known probe states |00>, |++>, |+i+i> are teleported through the noisy channel and the correction box, Bob
un-prepares them and measures: P(00) is the probe's fidelity. Under the twirl (``telecorr.twirl``) the probe along
axis a returns on each leg with probability f_q = (1 + lambda_q,a) / 2, lambda_q = diag(R_q T_q), so

    P(b5 b4) = product of the leg factors,        cost C = 1 - (1/3) sum_a P_a(00)     (Theorem T3.5).

The cost is pinned to the Haar-average fidelity (Theorems T3.6 / T5.9, the "sandwich"):

    (5/6)(1 - F_avg)  <=  C  <=  (5/6)(1 - F_avg) + (25/48)(1 - F_avg)^2,

so driving C to its minimum drives F_avg to its maximum, for every message, without knowing the message or the noise.
Shot noise: with N shots per probe the counts are multinomial; the measured cost is a sum of three binomials.
"""
from __future__ import annotations

from typing import Mapping, Sequence

import numpy as np

from .paulis import KEYS4
from .rotations import correction_rotations
from .twirl import leg_lambdas

#: probe name -> Bloch axis it tests (x = 0, y = 1, z = 2), in the order the circuits are run
PROBE_AXIS: dict[str, int] = {"|00>": 2, "|++>": 0, "|+i+i>": 1}
#: the quaternion form of the cost (Proposition P4.1): C = 1 - u^T M v / 3 with u, v the squared quaternions
M_QUAT = np.array([[3, 1, 1, 1], [1, 1, 0, 0], [1, 0, 1, 0], [1, 0, 0, 1]], float)


def probe_probs(params: Sequence[float] | None, T: Mapping[int, np.ndarray],
                probes: Sequence[str] = tuple(PROBE_AXIS)) -> np.ndarray:
    """Exact outcome probabilities (rows: probes; columns: keys "00", "01", "10", "11" = b5 b4)."""
    R = correction_rotations(params)
    lam = {q: leg_lambdas(R[q], T[q]) for q in (2, 3)}
    out = []
    for name in probes:
        a = PROBE_AXIS[name]
        f2, f3 = (1 + lam[2][a]) / 2, (1 + lam[3][a]) / 2            # leg q2 -> Bob q4 (bit b4), q3 -> q5 (b5)
        out.append([f3 * f2, f3 * (1 - f2), (1 - f3) * f2, (1 - f3) * (1 - f2)])
    return np.clip(np.array(out), 0, 1)


def cost(params: Sequence[float] | None, T: Mapping[int, np.ndarray], probes: Sequence[str] = tuple(PROBE_AXIS)) -> float:
    """Exact probe cost C = mean over the probes of 1 - P(00)."""
    return float(1 - probe_probs(params, T, probes)[:, 0].mean())


def sample_counts(rng: np.random.Generator, probs: np.ndarray, shots: int) -> np.ndarray:
    """Counts of each probe circuit (rows) drawn from the exact multinomial law."""
    return np.array([rng.multinomial(shots, p / p.sum()) for p in probs])


def counts_dicts(counts: np.ndarray) -> list[dict[str, int]]:
    return [{k: int(n) for k, n in zip(KEYS4, row) if n} for row in counts]


def sandwich(F_avg: float) -> tuple[float, float]:
    """Lower and upper bound of the probe cost for a given Haar-average fidelity (Theorem T5.9)."""
    g = 1 - F_avg
    return 5 / 6 * g, 5 / 6 * g + 25 / 48 * g * g


def quaternion_cost(u: np.ndarray, v: np.ndarray) -> float:
    """Probe cost from the squared quaternions u, v of the two residual rotations (Proposition P4.1)."""
    return float(1 - u @ M_QUAT @ v / 3)


def pauli_fingerprint(counts: np.ndarray, shots: int) -> dict[int, np.ndarray]:
    """Pauli rates (p_I, p_X, p_Y, p_Z) of each leg read from the three probe counts (probe order |00>, |++>, |+i+i>):
    f_a = p_I + p_a, so p_I = (f_X + f_Y + f_Z - 1) / 2 and p_a = f_a - p_I."""
    out = {}
    for q, bit in ((2, 1), (3, 0)):                      # key "b5b4": leg q2 is the last character
        f = {}
        for name, row in zip(PROBE_AXIS, counts):
            zero = sum(n for k, n in zip(KEYS4, row) if k[bit] == "0")
            f[PROBE_AXIS[name]] = zero / shots
        pI = (f[0] + f[1] + f[2] - 1) / 2
        out[q] = np.array([pI, f[0] - pI, f[1] - pI, f[2] - pI])
    return out
