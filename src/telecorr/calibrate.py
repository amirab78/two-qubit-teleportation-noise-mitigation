"""Calibrating the correction from probe counts only: COBYLA, the automatic recalibrator, and Rotosolve.

Roadmap: step 3 (COBYLA), step 4 (automatic recalibration), step 5 (Rotosolve in local coordinates).

The calibrator never sees the channel: it chooses six correction angles, a :class:`ProbeSampler` runs the three
probe circuits and returns counts, and the calibrator only reads the measured cost (``telecorr.probes``). The
sampler draws the counts from the exact law of the protocol (multinomial with the probabilities of T3.5), which is
the distribution a simulator or a device would produce; ``protocol.bob_counts`` gives the same distribution with Aer.

* :func:`cobyla` - COBYLA in phases of growing shots and shrinking steps (coarse, then fine; Stage 3).
* :class:`AutoCalibrator` - monitor -> detect (mean probe P(00) below a threshold) -> multi-start recalibration
  -> verify (Stage 4). An accepted start with C < 0.01 guarantees F > 0.985 for every message (P4.10).
* :func:`rotosolve` - every angle of a local coordinate enters the cost as an exact sinusoid; three costs fix it and
  the step jumps to its minimum (P5.15). The coordinates are local around the current correction,
  V = V_ref Rz(a) Ry(b) Rx(c) on each qubit, which are well conditioned everywhere.
* :class:`RelativeAutoCalibrator` - the step-5 calibrator for channels whose best P(00) is below 1: the trigger
  fires on a significant *change* of the monitored P(00), and the multi-start search stops when two starts agree.
"""
from __future__ import annotations

from typing import Mapping, Sequence

import numpy as np
from scipy.optimize import minimize

from .probes import PROBE_AXIS, pauli_fingerprint, probe_probs, sample_counts
from .rotations import angles_of, axis_rotation, correction_rotations, rotation_angle_axis, u_matrix
from .twirl import f_avg_of, message_fidelity

#: (shots per probe, initial step, max evaluations) of the COBYLA phases (Stages 3 and 4)
COBYLA_PHASES = ((2000, 0.6, 200), (20000, 0.15, 150), (100000, 0.04, 100))
COARSE, FINE = COBYLA_PHASES[:1], COBYLA_PHASES[1:]
#: (shots per probe, min sweeps, max sweeps) of the Rotosolve phases (Stage 5)
ROTO_COARSE, ROTO_FINE = ((2000, 2, 6),), ((20000, 2, 6), (100000, 1, 3))


class ProbeSampler:
    """The probe circuits of a channel with leg matrices T: exact-law counts, with a shot and evaluation budget.
    Set ``T`` to a new value to change the channel (drift, jumps); ``probes`` selects which probes are run."""

    def __init__(self, T: Mapping[int, np.ndarray], rng: np.random.Generator, probes: Sequence[str] = tuple(PROBE_AXIS)):
        self.T, self.rng, self.probes = dict(T), rng, tuple(probes)
        self.evals, self.shots = 0, 0

    def counts(self, params: Sequence[float] | None, shots: int) -> np.ndarray:
        probs = probe_probs(params, self.T, self.probes)
        self.evals += 1; self.shots += shots * len(probs)
        return sample_counts(self.rng, probs, shots)

    def cost(self, params: Sequence[float] | None, shots: int) -> float:
        """Measured cost: mean over the probes of 1 - n(00) / N."""
        return float(1 - self.counts(params, shots)[:, 0].mean() / shots)


# ---------------------------------------------------------------- true scores (for evaluation only)
def scores(params: Sequence[float] | None, T: Mapping[int, np.ndarray], messages: Sequence[np.ndarray]) -> dict:
    """Exact F of each message and the Haar average F_avg for a correction (the calibrator never uses this)."""
    R = correction_rotations(params)
    return {"F": [message_fidelity(psi, R, T) for psi in messages], "F_avg": f_avg_of(R, T)}


def residual_angles(params: Sequence[float] | None, T: Mapping[int, np.ndarray]) -> list[float]:
    """Rotation angle (rad) left on each leg by a correction, for a purely coherent channel (T a rotation)."""
    R = correction_rotations(params)
    return [rotation_angle_axis(R[q] @ T[q])[0] for q in (2, 3)]


# ---------------------------------------------------------------- COBYLA (Stages 3 and 4)
def cobyla(x0: Sequence[float], sampler: ProbeSampler, phases=COBYLA_PHASES, trace: list | None = None) -> np.ndarray:
    """Minimise the measured cost with COBYLA, phase by phase (more shots and smaller steps in each phase)."""
    x = np.array(x0, float)
    for shots, step, maxiter in phases:
        def f(p, s=shots):
            c = sampler.cost(p, s)
            if trace is not None:
                trace.append((sampler.shots, c, np.array(p)))
            return c
        x = minimize(f, x, method="COBYLA", options={"maxiter": maxiter, "rhobeg": step, "tol": 1e-5}).x
    return x


# ---------------------------------------------------------------- Rotosolve in local coordinates (Stage 5)
def recentre(x: Sequence[float], y: Sequence[float]) -> np.ndarray:
    """U(theta, phi, lambda) angles of V_ref(x) Rz(y0) Ry(y1) Rx(y2) on each channel qubit."""
    out = []
    for i in range(2):
        a, b, c = y[3 * i:3 * i + 3]
        V = u_matrix(*x[3 * i:3 * i + 3]) @ axis_rotation(2, a) @ axis_rotation(1, b) @ axis_rotation(0, c)
        out += list(angles_of(V))
    return np.array(out)


def rotosolve(x0: Sequence[float], sampler: ProbeSampler, phases=ROTO_COARSE + ROTO_FINE, trace: list | None = None) -> np.ndarray:
    """Coordinate descent with exact sinusoid steps; a phase stops when a whole sweep no longer lowers the measured
    cost by more than 2 sigma of its shot noise (after ``min_sweeps``)."""
    x = np.array(x0, float)
    n_probes = len(PROBE_AXIS)
    for shots, min_sweeps, max_sweeps in phases:
        cost = lambda y: sampler.cost(recentre(x, y), shots)                  # noqa: E731
        previous = None
        for sweep in range(max_sweeps):
            y = np.zeros(6)
            c0 = cost(y)
            sig = np.sqrt(2 * max(c0 * (1 - c0), 1e-3) / (n_probes * shots))
            if sweep >= min_sweeps and previous - c0 < max(2 * sig, 1e-4):
                break
            previous = c0
            for d in range(6):
                if d:
                    c0 = cost(y)
                yp = y.copy(); yp[d] += np.pi / 2
                ym = y.copy(); ym[d] -= np.pi / 2
                cp, cm = cost(yp), cost(ym)
                y[d] += np.arctan2(cm - cp, cp + cm - 2 * c0)               # minimum of a + b cos t + d sin t
                if trace is not None:
                    trace.append((sampler.shots, recentre(x, y)))
            x = recentre(x, y)
    return x


OPTIMISERS = {"cobyla": (cobyla, COARSE, FINE), "rotosolve": (rotosolve, ROTO_COARSE, ROTO_FINE)}


# ---------------------------------------------------------------- the automatic calibrator (Stage 4)
class AutoCalibrator:
    """Monitor the probe P(00); below ``threshold`` recalibrate by a multi-start search (previous angles, zero and
    random starts, each refined coarsely and judged with 20 000 shots per probe), stop early at an accepted start
    (cost < ``accept_cost``), refine the best one with the fine phases, and verify."""

    judge_shots = 20000

    def __init__(self, sampler: ProbeSampler, threshold: float = 0.99, monitor_shots: int = 4000,
                 n_random_starts: int = 6, accept_cost: float = 0.01, seed: int = 1, optimiser: str = "cobyla"):
        self.sampler, self.threshold, self.monitor_shots = sampler, threshold, monitor_shots
        self.n_random_starts, self.accept_cost = n_random_starts, accept_cost
        self.rng = np.random.default_rng(seed)
        self.optimiser = optimiser
        self.params = np.zeros(6)
        self.log: list[dict] = []

    def monitor(self) -> float:
        return 1 - self.sampler.cost(self.params, self.monitor_shots)

    def _starts(self) -> list[tuple[str, np.ndarray]]:
        starts = ([("previous angles", self.params.copy()), ("zero", np.zeros(6))] if np.any(self.params)
                  else [("zero (fresh calibrator)", np.zeros(6))])
        r = self.rng
        starts += [(f"random #{k + 1}", np.r_[r.uniform(0, np.pi), r.uniform(-np.pi, np.pi, 2),
                                              r.uniform(0, np.pi), r.uniform(-np.pi, np.pi, 2)])
                   for k in range(self.n_random_starts)]
        return starts

    def _stop(self, c: float, previous_best: float | None) -> str | None:
        return "cost below accept_cost" if c < self.accept_cost else None

    def recalibrate(self) -> dict:
        run, coarse, fine = OPTIMISERS[self.optimiser]
        e0, s0 = self.sampler.evals, self.sampler.shots
        tried, best, reason = [], None, "all starts tried"
        for label, x0 in self._starts():
            x = run(x0, self.sampler, coarse)
            c = self.sampler.cost(x, self.judge_shots)
            tried.append((label, c))
            previous_best = None if best is None else best[1]
            if best is None or c < best[1]:
                best = (x, c, label)
            why = self._stop(c, previous_best)
            if why:
                reason = why
                break
        self.params = run(best[0], self.sampler, fine)
        return {"starts": tried, "chosen": best[2], "reason": reason,
                "evaluations": self.sampler.evals - e0, "shots": self.sampler.shots - s0}

    def _fire(self, before: float) -> bool:
        return before < self.threshold

    def step(self) -> dict:
        """One monitoring cycle: measure; recalibrate and verify if the trigger fires."""
        before = self.monitor()
        entry = {"p00_before": before, "recalibrated": False}
        if self._fire(before):
            entry.update(recalibrated=True, **self.recalibrate())
            entry["p00_after"] = self.monitor()
        self.log.append(entry)
        return entry


class RelativeAutoCalibrator(AutoCalibrator):
    """The step-5 calibrator (Rotosolve by default). Trigger: the monitored P(00) differs from the reference verified
    after the last recalibration by more than max(min_margin, 3 sigma). Stop rule: an accepted start, a start that
    matches the floor of the last recalibration, or two starts that agree within ``agree_tol``. After recalibrating,
    a verification with ``verify_shots`` sets the new reference and reads the Pauli fingerprint of each leg."""

    def __init__(self, sampler: ProbeSampler, min_margin: float = 0.005, agree_tol: float = 0.01,
                 verify_shots: int = 20000, optimiser: str = "rotosolve", **kw):
        super().__init__(sampler, optimiser=optimiser, **kw)
        self.min_margin, self.agree_tol, self.verify_shots = min_margin, agree_tol, verify_shots
        self.reference: float | None = None
        self.floor: float | None = None
        self.fingerprint: dict | None = None

    def margin(self) -> float:
        r, n = self.reference, len(PROBE_AXIS)
        sig = np.sqrt(r * (1 - r) / (n * self.monitor_shots) + r * (1 - r) / (n * self.verify_shots))
        return max(self.min_margin, 3 * sig)

    def _stop(self, c, previous_best):
        if c < self.accept_cost:
            return "cost below accept_cost"
        if self.floor is not None and c <= self.floor + self.agree_tol:
            return "matches the floor of the last calibration"
        if previous_best is not None and abs(c - previous_best) <= self.agree_tol:
            return "two starts agree"
        return None

    def _fire(self, before):
        return self.reference is None or abs(before - self.reference) > self.margin()

    def recalibrate(self) -> dict:
        info = super().recalibrate()
        counts = self.sampler.counts(self.params, self.verify_shots)
        p00 = float(counts[:, 0].mean() / self.verify_shots)
        self.reference, self.floor = p00, 1 - p00
        self.fingerprint = pauli_fingerprint(counts, self.verify_shots)
        info["shots"] = info["shots"] + 3 * self.verify_shots
        info["p00_verified"] = p00
        return info
