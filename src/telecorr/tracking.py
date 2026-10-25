"""Following a drifting channel with few shots: a Kalman tracker on T and the ABC decision rule.

Roadmap: step 7 (efficient learning and tracking).

TOMO-9 learns T from scratch at every step. When the channel drifts slowly, yesterday's estimate is still worth a lot:
a Kalman filter keeps a Gaussian belief about vec(T) of each leg and updates it with every probe round (each round
measures three linear functionals of T; ``telecorr.tomography``). Drift enters as process noise in the rotation
directions of T, with a drift rate learned online from the innovations; a large normalised innovation (NIS above
the chi-square(6) 1e-4 point) declares a jump and restarts from a fresh TOMO-9 round.

ABC ("adaptive Bayesian calibration") then decides how much to measure: among 1, 2 or 3 rounds of the most
informative settings, or a fresh TOMO-9, each with N shots per circuit from a ladder, it takes the CHEAPEST option
whose predicted gap to the ceiling, (4/5) sum tr(H P)/2 over the legs, is below a target (Propositions P7.8, P7.9).
Information balance predicts the saving: gain = 1 + (target loss) / (one-step drift loss), large for slow drift.
"""
from __future__ import annotations

import numpy as np

from .ceiling import procrustes_batch
from .tomography import DESIGNS

CYC = DESIGNS["tomo9"]                                       # settings I, C, C^2
N_LADDER = np.array([100, 200, 500, 1000, 2000, 5000, 10000, 20000])
NIS_THRESHOLD = 27.9                                         # chi-square, 6 dof (2 legs x 3 circuits), p = 1e-4
SCENARIOS = {                                                # drift per step: angle sd (rad), axis sd; jumps; changes of p
    "very slow": dict(sig_ang=0.001, sig_ax=0.0015, jumps=(), p_changes=()),
    "slow": dict(sig_ang=0.005, sig_ax=0.0075, jumps=(), p_changes=()),
    "medium": dict(sig_ang=0.015, sig_ax=0.0225, jumps=(20, 40), p_changes=(10, 30, 50)),
    "fast": dict(sig_ang=0.04, sig_ax=0.06, jumps=(20, 40), p_changes=(10, 30, 50)),
}


# ---------------------------------------------------------------- batched helpers
def trace_RT(R: np.ndarray, T: np.ndarray) -> np.ndarray:
    return np.einsum("nab,nba->n", R, T)


def f_avg2(tr2: np.ndarray, tr3: np.ndarray) -> np.ndarray:
    """F_avg of a channel from tr(R T) of its two legs."""
    return (4 * (1 + tr2) / 4 * (1 + tr3) / 4 + 1) / 5


def _skew(v: np.ndarray) -> np.ndarray:
    z = np.zeros(v.shape[:-1])
    return np.stack([np.stack([z, -v[..., 2], v[..., 1]], -1), np.stack([v[..., 2], z, -v[..., 0]], -1),
                     np.stack([-v[..., 1], v[..., 0], z], -1)], -2)


def rodrigues(v: np.ndarray) -> np.ndarray:
    """Rotation matrices of rotation vectors v (n, 3)."""
    ang = np.linalg.norm(v, axis=1)
    K = _skew(v / np.maximum(ang, 1e-12)[:, None]); a = ang[:, None, None]
    return np.eye(3) + np.sin(a) * K + (1 - np.cos(a)) * K @ K


def incoherent_diag_batch(fam: np.ndarray, s: np.ndarray) -> np.ndarray:
    """diag of the incoherent T for family index (0 none, 1 D, 2 A, 3 M) and strength s."""
    g = np.where(fam == 2, s, np.where(fam == 3, s / 2, 0.0))
    lam = np.where(fam == 3, s, 0.0)
    d = np.where(fam == 1, s, 0.0)
    a = np.sqrt(1 - g) * np.sqrt(1 - lam) * (1 - d)
    return np.stack([a, a, (1 - g) * (1 - d)], -1)


def measure_lambda(rng: np.random.Generator, R: np.ndarray, T: np.ndarray, N: np.ndarray) -> np.ndarray:
    """One probe round at setting R (n, 3, 3) on legs T (n, 3, 3): lambda_hat (n, 3) with binomial shot noise."""
    lam = np.einsum("nak,nka->na", R, T)
    return 2 * rng.binomial(N[:, None], np.clip((1 + lam) / 2, 0, 1)) / N[:, None] - 1


def round_G(S: np.ndarray) -> np.ndarray:
    """Linear functionals (..., 3, 9) of vec(T) measured by the 3 probe circuits of settings S (..., 3, 3)."""
    G = np.zeros(S.shape[:-2] + (3, 9))
    for a in range(3):
        G[..., a, [3 * c + a for c in range(3)]] = S[..., a, :]
    return G


A9 = round_G(CYC).reshape(9, 9)                              # the 9 circuits of TOMO-9: a permutation matrix
A9_inv = np.linalg.inv(A9)


def loss_hessian(T: np.ndarray, eps: float = 1e-3) -> np.ndarray:
    """Hessian (n, 9, 9) of the leg loss L(T_hat) = (tr(R*(T) T) - tr(procrustes(T_hat) T)) / 4 at T_hat = T."""
    n = len(T); best = trace_RT(procrustes_batch(T), T)
    L = lambda E: (best - trace_RT(procrustes_batch(T + E.reshape(1, 3, 3)), T)) / 4     # noqa: E731
    H = np.zeros((n, 9, 9))
    for i in range(9):
        for j in range(i, 9):
            ei = np.zeros(9); ei[i] = eps; ej = np.zeros(9); ej[j] = eps
            H[:, i, j] = H[:, j, i] = (L(ei + ej) - L(ei - ej) - L(ej - ei) + L(-ei - ej)) / (4 * eps * eps)
    return H


# ---------------------------------------------------------------- the Kalman tracker
class KalmanTracker:
    """Kalman filters on vec(T) for M legs (legs 2i, 2i + 1 = q2, q3 of channel i). Process noise in the rotation
    directions of T on both sides, Q = s_rot sum_k [vec(T E_k) vec(T E_k)^T + vec(E_k T) vec(E_k T)^T] + q_iso I,
    with the drift rate s_rot learned from the innovations (multiplicative steps)."""
    E = np.array([[[0, 0, 0], [0, 0, -1], [0, 1, 0]], [[0, 0, 1], [0, 0, 0], [-1, 0, 0]],
                  [[0, -1, 0], [1, 0, 0], [0, 0, 0]]], float)

    def __init__(self, M: int, s0: float = 1e-4, q_iso: float = 1e-7):
        self.M, self.q_iso = M, q_iso
        self.x = np.zeros((M, 9)); self.P = np.broadcast_to(np.eye(9), (M, 9, 9)).copy()
        self.s = np.full(M, s0)

    def T_hat(self) -> np.ndarray:
        return self.x.reshape(self.M, 3, 3)

    def _Qrot(self) -> np.ndarray:
        T = self.T_hat()
        V = np.concatenate([np.einsum("mab,kbc->mkac", T, self.E), np.einsum("kab,mbc->mkac", self.E, T)], 1).reshape(self.M, 6, 9)
        return np.einsum("mki,mkj->mij", V, V)

    def init_tomo9(self, lam: np.ndarray, N: np.ndarray, idx: np.ndarray | None = None) -> None:
        """(Re)start legs idx from one TOMO-9 round: least-squares mean and its covariance."""
        idx = np.arange(self.M) if idx is None else idx
        x = lam @ A9_inv.T
        v = np.maximum(1 - (x @ A9.T) ** 2, 1.0 / N[:, None]) / N[:, None]
        self.x[idx] = x
        self.P[idx] = np.einsum("ij,mj,kj->mik", A9_inv, v, A9_inv)

    def predict(self) -> None:
        self._Q = self._Qrot()
        self.P = self.P + self.s[:, None, None] * self._Q + self.q_iso * np.eye(9)

    def update(self, S: np.ndarray, lam: np.ndarray, N: np.ndarray, active: np.ndarray | None = None) -> np.ndarray:
        """Sequential scalar updates with one round (settings S, results lam, N shots); returns the NIS per leg."""
        active = np.ones(self.M, bool) if active is None else active
        G = round_G(S); nis = np.zeros(self.M); Nn = np.maximum(N, 1)
        for a in range(3):
            g = G[:, a]
            pred = np.einsum("mi,mi->m", g, self.x)
            v = np.maximum(1 - np.clip(pred, -1, 1) ** 2, 1.0 / Nn) / Nn
            Pg = np.einsum("mij,mj->mi", self.P, g)
            s = np.einsum("mi,mi->m", g, Pg) + v
            nu = np.where(active, lam[:, a] - pred, 0.0)
            K = np.where(active[:, None], Pg / s[:, None], 0.0)
            self.x = self.x + K * nu[:, None]
            self.P = self.P - K[:, :, None] * Pg[:, None, :]
            nis += np.where(active, nu ** 2 / s, 0.0)
        self.P = 0.5 * (self.P + np.swapaxes(self.P, 1, 2))
        return nis

    def adapt(self, mask: np.ndarray, nis1: np.ndarray, eta: float = 0.2) -> None:
        """Innovation matching: the NIS of a round has mean 3 if the drift model is right."""
        step = np.clip(eta * (nis1 / 3.0 - 1.0), -0.5, 0.5)
        self.s = np.where(mask, np.clip(self.s * np.exp(step), 1e-9, 1.0), self.s)


def _post_cov(P, G, x, N):
    lam = np.einsum("mai,mi->ma", G, x)
    v = np.maximum(1 - np.clip(lam, -1, 1) ** 2, 1.0 / N[:, None]) / N[:, None]
    return np.linalg.inv(np.linalg.inv(P + 1e-12 * np.eye(9)) + np.einsum("mai,ma,maj->mij", G, 1 / v, G))


def predicted_gap(P: np.ndarray, H: np.ndarray, pI: np.ndarray) -> np.ndarray:
    """Predicted two-leg gap per channel: (4/5)(p_I(q3) L2 + p_I(q2) L3), L = tr(H P)/2 per leg."""
    L = 0.5 * np.einsum("mij,mji->m", H, P)
    return (4 / 5) * (pI[1::2] * L[0::2] + pI[0::2] * L[1::2])


def abc_choose(kf: KalmanTracker, H: np.ndarray, pI: np.ndarray, target: float, N_nominal: int = 1000, max_rounds: int = 3):
    """The ABC rule: the cheapest option (1-3 greedy rounds, or a fresh TOMO-9; N from the ladder) whose predicted
    gap is below ``target``; if none is, the lowest predicted gap at the largest N.
    Returns settings (M, max_rounds, 3, 3), N (M,), rounds (M,), fresh (M,)."""
    M, n = kf.M, kf.M // 2
    Rh = procrustes_batch(kf.T_hat())
    cands = np.concatenate([np.broadcast_to(CYC, (M, 3, 3, 3)), np.einsum("kab,mbc->mkac", CYC, Rh)], 1)
    chosen = np.zeros((M, max_rounds, 3, 3)); P_cur = kf.P.copy()
    for r in range(max_rounds):                                # greedy settings, judged at the nominal N
        best, bestS, bestP = np.full(M, np.inf), np.zeros((M, 3, 3)), P_cur.copy()
        for k in range(cands.shape[1]):
            Pk = _post_cov(P_cur, round_G(cands[:, k]), kf.x, np.full(M, N_nominal))
            Lk = 0.5 * np.einsum("mij,mji->m", H, Pk)
            better = Lk < best
            best = np.where(better, Lk, best); bestS[better] = cands[better, k]; bestP[better] = Pk[better]
        chosen[:, r] = bestS; P_cur = bestP
    cost = np.full(n, np.inf); pick_N = np.full(n, N_LADDER[-1]); pick_r = np.full(n, max_rounds)
    fallback = np.full(n, np.inf); fresh = np.zeros(n, bool)
    G9 = round_G(np.broadcast_to(CYC, (M, 3, 3, 3))).reshape(M, 9, 9)
    lam9 = np.einsum("mci,mi->mc", G9, kf.x)
    for N in N_LADDER:                                         # option "fresh TOMO-9": forget the memory
        v = np.maximum(1 - np.clip(lam9, -1, 1) ** 2, 1.0 / N) / N
        pg = predicted_gap(np.linalg.inv(np.einsum("mci,mc,mcj->mij", G9, 1 / v, G9)), H, pI)
        ok = (pg <= target) & (9 * N < cost)
        cost[ok] = 9 * N; pick_N[ok] = N; pick_r[ok] = 3; fresh[ok] = True
        if N == N_LADDER[-1]:
            fallback = np.where(~np.isfinite(cost), pg, fallback); fresh |= ~np.isfinite(cost)
    for r in range(1, max_rounds + 1):
        for N in N_LADDER:
            Pn = kf.P
            for j in range(r):
                Pn = _post_cov(Pn, round_G(chosen[:, j]), kf.x, np.full(M, N))
            pg = predicted_gap(Pn, H, pI)
            ok = (pg <= target) & (3 * r * N < cost)
            cost[ok] = 3 * r * N; pick_N[ok] = N; pick_r[ok] = r; fresh[ok] = False
            if N == N_LADDER[-1]:
                worse = ~np.isfinite(cost) & (pg < fallback)
                fallback[worse] = pg[worse]; pick_r[worse] = r; pick_N[worse] = N; fresh[worse] = False
    fr = np.repeat(fresh, 2)
    chosen[fr] = np.broadcast_to(CYC, (int(fr.sum()), 3, 3, 3))
    return chosen, np.repeat(pick_N, 2), np.repeat(pick_r, 2), fr


def run_abc(Ts: np.ndarray, target: float, seed: int = 0, N_init: int = 3000, N_reset: int = 3000) -> dict:
    """ABC over channel sequences Ts (steps, n, 2, 3, 3): TOMO-9 start, then predict -> choose -> measure -> update
    -> jump test -> correction = Procrustes(estimate). Returns the exact gap and the shots of every step."""
    rng = np.random.default_rng(seed)
    steps, n = Ts.shape[:2]; M = 2 * n
    kf = KalmanTracker(M)
    gaps, shots, jumps = np.zeros((steps, n)), np.zeros((steps, n)), np.zeros((steps, n), bool)
    for t in range(steps):
        T = Ts[t].reshape(M, 3, 3)
        if t == 0:
            Nv = np.full(M, N_init)
            lam = np.stack([measure_lambda(rng, np.broadcast_to(S, (M, 3, 3)), T, Nv) for S in CYC], 1).reshape(M, 9)
            kf.init_tomo9(lam, Nv); shots[t] = 9 * N_init
        else:
            kf.predict()
            That = kf.T_hat()
            H = loss_hessian(That)
            pI = (1 + trace_RT(procrustes_batch(That), That)) / 4
            S, Nv, R_n, fr = abc_choose(kf, H, pI, target)
            kf.P[fr] = 10.0 * np.eye(9)                         # "fresh": diffuse prior, then TOMO-9
            nis1 = None
            for r in range(S.shape[1]):
                active = R_n > r
                if not active.any():
                    break
                nis = kf.update(S[:, r], measure_lambda(rng, S[:, r], T, np.maximum(Nv, 1)), Nv, active)
                nis1 = nis if r == 0 else nis1
            jump = (nis1[0::2] + nis1[1::2] > NIS_THRESHOLD) & ~fr[0::2]
            kf.adapt(~np.repeat(jump, 2) & ~fr, nis1)
            shots[t] = 3 * Nv[0::2] * R_n[0::2]
            if jump.any():
                idx = np.where(np.repeat(jump, 2))[0]
                Nr = np.full(len(idx), N_reset)
                lamr = np.stack([measure_lambda(rng, np.broadcast_to(S_, (len(idx), 3, 3)), T[idx], Nr) for S_ in CYC], 1)
                kf.init_tomo9(lamr.reshape(len(idx), 9), Nr, idx)
                shots[t, jump] += 9 * N_reset; jumps[t] = jump
        R = procrustes_batch(kf.T_hat())
        b = trace_RT(procrustes_batch(T), T)
        gaps[t] = f_avg2(b[0::2], b[1::2]) - f_avg2(trace_RT(R, T)[0::2], trace_RT(R, T)[1::2])
    return dict(gap=gaps, shots=shots, jumps=jumps)


def run_tomo9(Ts: np.ndarray, N: int, seed: int = 0) -> dict:
    """The baseline: a fresh TOMO-9 (9 N shots) at every step."""
    from .tomography import ls_T
    rng = np.random.default_rng(seed)
    steps, n = Ts.shape[:2]; M = 2 * n
    gaps = np.zeros((steps, n))
    S = np.broadcast_to(CYC, (M, 3, 3, 3))
    for t in range(steps):
        T = Ts[t].reshape(M, 3, 3); Nv = np.full(M, N)
        lam = np.stack([measure_lambda(rng, S[:, j], T, Nv) for j in range(3)], 1)
        R = procrustes_batch(ls_T(S, lam))
        b = trace_RT(procrustes_batch(T), T)
        gaps[t] = f_avg2(b[0::2], b[1::2]) - f_avg2(trace_RT(R, T)[0::2], trace_RT(R, T)[1::2])
    return dict(gap=gaps, shots=np.full((steps, n), 9 * N))


# ---------------------------------------------------------------- drifting channels
def drift_scenarios(rng: np.random.Generator, n: int, steps: int = 60, sig_ang: float = 0.04, sig_ax: float = 0.06,
                    jumps=(20, 40), p_changes=(10, 30, 50)) -> np.ndarray:
    """n episodes of two-leg channels (steps, n, 2, 3, 3): rotation axes and angles random-walk, new random rotations
    at the jumps, family D / A / M at random with p ~ U(0.02, 0.2) redrawn at p_changes (q3 at 1.5 p)."""
    out = np.zeros((steps, n, 2, 3, 3))
    axis = rng.normal(size=(n, 2, 3)); angle = rng.uniform(0, np.pi, (n, 2))
    fam = rng.integers(1, 4, n); p = rng.uniform(0.02, 0.2, n)
    rot = lambda ax, an: rodrigues((ax / np.linalg.norm(ax, axis=-1, keepdims=True) * an[..., None]).reshape(-1, 3)).reshape(n, 2, 3, 3)  # noqa: E731
    for t in range(steps):
        if t in jumps:
            axis = rng.normal(size=(n, 2, 3)); angle = rng.uniform(0, np.pi, (n, 2))
        elif t > 0:
            axis = axis / np.linalg.norm(axis, axis=-1, keepdims=True) + rng.normal(0, sig_ax, (n, 2, 3))
            angle = np.clip(angle + rng.normal(0, sig_ang, (n, 2)), 0, np.pi)
        if t in p_changes:
            p = rng.uniform(0.02, 0.2, n)
        R = rot(axis, angle)
        out[t, :, 0] = incoherent_diag_batch(fam, p)[:, :, None] * R[:, 0]
        out[t, :, 1] = incoherent_diag_batch(fam, 1.5 * p)[:, :, None] * R[:, 1]
    return out


def summary(res: dict) -> dict:
    g, s = res["gap"], res["shots"]
    return dict(mean_gap=float(g.mean()), p99_gap=float(np.percentile(g, 99)), shots_per_step=float(s[1:].mean()),
                constant=float(g.mean() * s[1:].mean()))
