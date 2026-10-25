"""The hardware analysis of the two rounds on ``ibm_kingston`` (4 October 2026), from the raw counts.

Roadmap: step 8 (real hardware).

Round 1 measured the device (readout calibration, native teleportation in versions A / B / C, the Alice reference),
the injected Stage 2 noise N1 without correction and with the exact inverse ("oracle"), an 8 us wait (N2), and
the TOMO-9 / TOMO-12 probe circuits that learn the correction from counts. Round 2 ran the learned corrections.

Block names: ``R1_nativeA`` (version A, Pauli frame), ``R1_nativeB`` (B, mid-circuit measurement + feed-forward),
``R1_C`` (C, deferred measurement), ``R1_N1none`` / ``R1_N1oracle`` / ``R1_N2none``, ``R2_N1corr`` (N1 + TOMO-9 from
real counts), ``R2_N2corr``, ``R2_N0corr`` (native + TOMO-9 of the device fingerprint), ``R2_DDN2none`` /
``R2_DDN2corr`` (with dynamical decoupling); ``...B`` = version B.
"""
from __future__ import annotations

import datetime as dt
import math
from typing import Mapping, Sequence

import numpy as np

from .. import tomography as proc
from ..ceiling import procrustes
from ..detector import detect
from ..paulis import TOMO_BASES
from ..rotations import oracle_angles, so3, u_matrix
from ..states import MESSAGES, message_state
from ..twirl import leg_state_fidelity
from . import budget as B
from . import data as D
from .frame import frame_counts, frame_keys, raw_counts
from .readout import READOUT_PATTERNS, confusion_from_counts, mitigate_bob
from .tomography import analyse, linear_F, linear_F_variance, linear_weights

#: the configuration fixed before Round 1 (8c); the values in 8c_choices.json override the defaults
DEFAULTS = {"B": "B", "tau_ns": 6000, "shots_tomo": 3000, "shots_cal": 2000, "C": True,
            "dd_blocks": ["N2none", "N2corr"], "twirl_blocks": [], "n0corr_rule_rad": 0.05}


def config() -> dict:
    return dict(DEFAULTS, **D.preregistration("8c_choices.json")["cfg"])


# --------------------------------------------------------------------------------------- count tables
def counts9(data: Mapping, prefix: str, state: str, version: str = "A", M=None) -> dict:
    """Bob's 9 tomography tables of one block (Pauli frame for A; optional Bob-only readout mitigation)."""
    out = {}
    for bs in TOMO_BASES:
        d = data[f"{prefix}_{state}_{bs}"]["ints"]
        a, b = np.array(d["alice"]), np.array(d["bob"])
        if M is not None:
            q = mitigate_bob(a, b, bs, M, version)
            out[bs] = {k: v * len(b) for k, v in q.items()}
        elif version == "A":
            out[bs] = frame_counts(a, b, bs)
        else:
            out[bs] = raw_counts(b)
    return out


def probe_dists(data: Mapping, labels: Sequence[str]) -> list[dict]:
    """Frame-corrected counts of TOMO probe circuits; each label ends with the probe axis."""
    out = []
    for lab in labels:
        d = data[lab]["ints"]
        out.append(frame_counts(np.array(d["alice"]), np.array(d["bob"]), lab[-1] * 2))
    return out


def unprep_P00(data: Mapping, label: str) -> tuple[float, float]:
    """P(00) of an un-preparation circuit (version B) and its binomial sigma."""
    b = np.array(data[label]["ints"]["bob"])
    p = float(np.mean(b == 0))
    return p, float(np.sqrt(p * (1 - p) / len(b)))


# --------------------------------------------------------------------------------------- the two rounds
ROUND1_BLOCKS = [("R1_nativeA", "A", MESSAGES), ("R1_nativeB", "B", MESSAGES), ("R1_N1none", "A", MESSAGES),
                 ("R1_N1oracle", "A", MESSAGES), ("R1_N2none", "A", MESSAGES), ("R1_N1noneB", "B", ("bell_30_70",))]


def analyse_round1(data: Mapping, cfg: dict | None = None, boot: int = 300, mitigate: bool = True) -> dict:
    """Round 1: readout calibration, every tomography block (raw and Bob-mitigated), un-preparation checks, TOMO-9 /
    TOMO-12 of the device fingerprint, of N1 and of N2, the crosstalk detector, and the Alice reference."""
    cfg = cfg or config()
    res: dict = {"tomo": {}, "unprep": {}, "tomo_mit": {}}
    res["readout"] = confusion_from_counts({pat: data[f"readout_{pat}"]["counts"] for pat in READOUT_PATTERNS})
    M = res["readout"]["M_bob"] if mitigate else None
    blocks = list(ROUND1_BLOCKS) + ([("R1_C", "B", ("product_HT",))] if cfg["C"] else [])
    for prefix, ver, states in blocks:
        for st in states:
            psi = message_state(st)
            res["tomo"][f"{prefix}_{st}"] = analyse(counts9(data, prefix, st, ver), psi, boot)
            if M is not None:
                res["tomo_mit"][f"{prefix}_{st}"] = analyse(counts9(data, prefix, st, ver, M), psi, 0)
    for lab in ("R1_nativeB_unprep", "R1_N1none_unprep", "R1_N1oracle_unprep"):
        res["unprep"][lab] = unprep_P00(data, lab)
    lab9, lab12 = proc.probe_labels("tomo9"), proc.probe_labels("tomo12")
    res["fingerprint"] = proc.estimate(probe_dists(data, ["R1_fp_" + l for l in lab9]), "tomo9")
    res["cal9"] = proc.estimate(probe_dists(data, ["R1_cal9_" + l for l in lab9]), "tomo9")
    res["cal12"] = proc.estimate(probe_dists(data, ["R1_cal12_" + l for l in lab12]), "tomo12")
    res["calN2"] = proc.estimate(probe_dists(data, ["R1_calN2_" + l for l in lab9]), "tomo9")
    res["detector_fp"] = detect(probe_dists(data, ["R1_fp_" + l for l in lab9]), cfg["shots_cal"])
    res["detector_cal"] = detect(probe_dists(data, ["R1_cal9_" + l for l in lab9]), cfg["shots_cal"])
    res["alice"] = {}
    for st in MESSAGES:
        c9 = {bs: {k: int(v) for k, v in data[f"alice_{st}_{bs}"]["counts"]["msg"].items()} for bs in TOMO_BASES}
        res["alice"][st] = {"counts9": c9, **analyse(c9, message_state(st), 0)}
    return res


def corrections_from_round1(r1: Mapping, cfg: dict | None = None) -> dict:
    """The Round 2 corrections, with the pre-registered rule for the N0 block: N1 <- TOMO-9 of the injected noise,
    N2 <- TOMO-9 of the wait, N0 <- the device fingerprint only if its largest rotation exceeds 0.05 rad."""
    cfg = cfg or config()
    fp = r1["fingerprint"]
    ang = max(fp["rotation_2"]["angle"], fp["rotation_3"]["angle"])
    return {"N1": list(r1["cal9"]["x_new"]), "N2": list(r1["calN2"]["x_new"]),
            "N0": list(fp["x_new"]) if ang > cfg["n0corr_rule_rad"] else None, "fingerprint_angle": ang}


def analyse_round2(data: Mapping, corr: Mapping, M, cfg: dict | None = None, boot: int = 300) -> dict:
    """Round 2: every corrected block (raw and mitigated with Round 1's readout matrix M), the un-preparation check
    and the verification TOMO-9 measured after the N1 correction."""
    cfg = cfg or config()
    res: dict = {"tomo": {}, "unprep": {}, "tomo_mit": {}}
    blocks = [("R2_N1corr", "A", MESSAGES), ("R2_N2corr", "A", MESSAGES), ("R2_N1corrB", "B", ("bell_30_70",))]
    if corr.get("N0") is not None:
        blocks.append(("R2_N0corr", "A", MESSAGES))
    blocks += [(f"R2_DD{b}", "A", ("product_HT",)) for b in cfg["dd_blocks"]]
    blocks += [(f"R2_TW{b}", "A", ("product_HT",)) for b in cfg["twirl_blocks"]]
    for prefix, ver, states in blocks:
        for st in states:
            psi = message_state(st)
            res["tomo"][f"{prefix}_{st}"] = analyse(counts9(data, prefix, st, ver), psi, boot)
            if M is not None:
                res["tomo_mit"][f"{prefix}_{st}"] = analyse(counts9(data, prefix, st, ver, M), psi, 0)
    res["unprep"]["R2_N1corr_unprep"] = unprep_P00(data, "R2_N1corr_unprep")
    res["verification"] = proc.estimate(probe_dists(data, ["R2_ver_" + l for l in proc.probe_labels("tomo9")]),
                                        "tomo9", corr["N1"])
    return res


# --------------------------------------------------------------------------------------- the twirl model
MODEL_CASES = {"R1_nativeA": ("fingerprint", None), "R1_N1none": ("cal9", None), "R1_N1oracle": ("cal9", "oracle"),
               "R1_N2none": ("calN2", None), "R2_N0corr": ("fingerprint", "learned"), "R2_N1corr": ("cal9", "learned"),
               "R2_N2corr": ("calN2", "learned")}


def twirl_model(r1: Mapping) -> dict:
    """Prediction of every block from the Round 1 TOMO-9 matrices only: each leg's Bob channel is the twirled map
    diag(diag(R T)) with T = the TOMO-9 matrix of the block's channel (fingerprint for native, cal9 for N1, calN2 for
    N2) and R = no correction, the oracle, or the learned Procrustes rotation. Readout errors are inside T."""
    T = {k: B.leg_T(r1, k) for k in ("fingerprint", "cal9", "calN2")}
    Ropt = {k: {q: procrustes(T[k][q]) for q in (2, 3)} for k in T}
    ov = oracle_angles("N1")
    orc = {2: so3(u_matrix(*ov[:3])), 3: so3(u_matrix(*ov[3:]))}
    out: dict = {"blocks": {}}
    for key, (src, corr) in MODEL_CASES.items():
        for st in MESSAGES:
            Mq = {}
            for q in (2, 3):
                Rq = np.eye(3) if corr is None else (orc[q] if corr == "oracle" else Ropt[src][q])
                Mq[q] = np.diag(np.diag(Rq @ T[src][q]))
            F = leg_state_fidelity(message_state(st), Mq[2], Mq[3])
            hwF = r1["tomo"].get(f"{key}_{st}", {}).get("F") if key.startswith("R1_") else None
            out["blocks"][f"{key}_{st}"] = {"model_F": F, "hardware_F": hwF, "source_T": src, "correction": corr or "none"}
    chk = [v["model_F"] - v["hardware_F"] for v in out["blocks"].values() if v["hardware_F"] is not None]
    out["check_round1"] = {"mean_diff": float(np.mean(chk)), "max_abs_diff": float(np.max(np.abs(chk))), "n": len(chk)}
    return out


# --------------------------------------------------------------------------------------- per-circuit check
def circuit_marginals(entry: Mapping, label: str, version: str) -> dict[str, np.ndarray]:
    """Alice's 16-outcome and Bob's 4-outcome distributions of one circuit (Bob frame-corrected in version A)."""
    ints = entry["ints"]
    if version == "msg":
        m = np.asarray(ints["msg"])
        return {"msg": np.bincount(m, minlength=4) / len(m)}
    a = np.asarray(ints["alice"]); b = np.asarray(ints["bob"])
    out = {"alice": np.bincount(a, minlength=16) / len(a)}
    if version == "A":
        tail = label.split("_")[-1]
        bases = tail if (len(tail) == 2 and set(tail) <= set("XYZ")) else (tail * 2 if tail in ("X", "Y", "Z") else "ZZ")
        out["bob"] = np.bincount(frame_keys(a, b, bases), minlength=4) / len(b)
    else:
        out["bob"] = np.bincount(b, minlength=4) / len(b)
    return out


def tvd(p, q) -> float:
    return 0.5 * float(np.abs(np.asarray(p) - np.asarray(q)).sum())


def circuit_check(data: Mapping, pred: Mapping, band: Mapping, shots: Mapping, versions: Mapping, seed: int = 11,
                  margin: float = 0.05, reps: int = 200) -> list[dict]:
    """Every circuit's measured marginals against the pre-registered prediction: flagged when, for any register,
    TVD(measured, predicted) > (99.9 % quantile of the shot-noise TVD at its shots) + max(margin, model band)."""
    rng = np.random.default_rng(seed)
    rows = []
    for lab, entry in data.items():
        if lab not in pred:
            rows.append({"label": lab, "ok": False, "detail": "no prediction"}); continue
        m = circuit_marginals(entry, lab, versions[lab])
        n = shots.get(lab) or len(next(iter(entry["ints"].values())))
        worst, ok = None, True
        for reg, p in m.items():
            pp = np.asarray(pred[lab][reg])
            pc = np.clip(pp.astype(float), 0, None); pc = pc / pc.sum()
            thr0 = float(np.quantile([tvd(rng.multinomial(n, pc) / n, pc) for _ in range(reps)], 0.999))
            thr = thr0 + max(margin, band.get(lab, {}).get(reg, 0.0))
            d = tvd(p, pp)
            if worst is None or d - thr > worst[1] - worst[2]:
                worst = (reg, d, thr)
            ok &= d <= thr
        rows.append({"label": lab, "ok": bool(ok), "register": worst[0], "tvd": worst[1], "threshold": worst[2]})
    return rows


# --------------------------------------------------------------------------------------- hardware vs predictions
BLOCK_TEXT = {"R1_nativeA": "native, A (frame)", "R1_nativeB": "native, B (dynamic)", "R1_C": "native, C (deferred)",
              "R1_N1none": "N1, no correction (A)", "R1_N1noneB": "N1, no correction (B)", "R1_N1oracle": "N1 + oracle (A)",
              "R2_N1corr": "N1 + TOMO-9 from real counts (A)", "R2_N1corrB": "N1 + TOMO-9 from real counts (B)",
              "R1_N2none": "N2 (8 µs wait), no correction", "R2_N2corr": "N2 + TOMO-9", "R2_DDN2none": "N2 + DD",
              "R2_DDN2corr": "N2 + TOMO-9 + DD", "R2_N0corr": "native + TOMO-9 (N0 block)"}


# Plain names for figures and tables (version A unless stated). ``BLOCK_TEXT`` above is the wording of the saved results.
NAMES = {"R1_nativeA": "device alone", "R1_nativeB": "device alone, version B", "R1_C": "device alone, version C",
         "R1_N1none": "N1, no correction", "R1_N1noneB": "N1, no correction, version B",
         "R1_N1oracle": "N1, exact inverse (oracle)", "R2_N1corr": "N1, TOMO-9 correction",
         "R2_N1corrB": "N1, TOMO-9 correction, version B", "R2_N0corr": "device alone, TOMO-9 correction",
         "R1_N2none": "8 µs wait, no correction", "R2_N2corr": "8 µs wait, TOMO-9 correction",
         "R2_DDN2none": "8 µs wait, DD", "R2_DDN2corr": "8 µs wait, TOMO-9 correction + DD"}


def describe(key: str) -> tuple[str, int, str]:
    """``"R2_N1corr_bell_30_70"`` -> ("N1, TOMO-9 correction", 2, "bell_30_70"): plain name, round, message."""
    pre, st = _split_key(key)
    return NAMES.get(pre, pre), int(pre[1]), st


def _split_key(key: str) -> tuple[str, str]:
    st = "bell_30_70" if key.endswith("bell_30_70") else "product_HT"
    return key[: -len(st) - 1], st


def compare(r1: Mapping, r2: Mapping, P1: Mapping, P2: Mapping, A: Mapping | None) -> list[dict]:
    """Every hardware fidelity next to its time-stamped prediction (8c, and the Round 2 addendum): difference,
    z-score with the hardware and twin spreads, inside / outside the pre-registered band."""
    rows = []
    for rnd, R, P in ((1, r1, P1), (2, r2, P2)):
        for key, v in R["tomo"].items():
            pre, st = _split_key(key)
            pb = (P["blocks"].get(key) or {}).get("F") or {}
            row = {"round": rnd, "key": key, "block": BLOCK_TEXT.get(pre, pre), "state": st, "F": v["F"],
                   "F_sigma": v.get("F_sigma"), "F_mle": v.get("F_mle"), "purity": v.get("purity"),
                   "F_mit": (R.get("tomo_mit", {}).get(key) or {}).get("F"), "pred": pb.get("prediction"),
                   "pred_L3plus": (pb.get("L3plus") or {}).get("mean"), "band": pb.get("band"),
                   "pred_sd": (pb.get("L3") or {}).get("sd"),
                   "pred_mit": ((P["blocks"].get(key) or {}).get("F_mitigated") or {}).get("mean")}
            if A and rnd == 2 and key in A["blocks"]:
                ab = A["blocks"][key]["F"]
                row.update({"addendum": ab["prediction"], "addendum_L3plus": (ab.get("L3plus") or {}).get("mean"),
                            "addendum_band": ab["band"]})
            if row["pred"] is not None:
                s = math.sqrt((row["F_sigma"] or 0) ** 2 + (row["pred_sd"] or 0) ** 2) or 1e-9
                row["diff"] = row["F"] - row["pred"]; row["z"] = row["diff"] / s
                row["inside_band"] = bool(row["band"][0] <= row["F"] <= row["band"][1])
            rows.append(row)
    return rows


def hypotheses(r1: Mapping, r2: Mapping, P1: Mapping) -> list[dict]:
    """The pre-registered hypotheses H1-H11 of 8c evaluated on the real numbers (verdict wording fixed in 8c)."""
    F = lambda R, k: R["tomo"][k]["F"]                                         # noqa: E731
    S = lambda R, k: R["tomo"][k].get("F_sigma") or 0.0                        # noqa: E731
    H = []
    fp = max(r1["fingerprint"]["rotation_2"]["angle"], r1["fingerprint"]["rotation_3"]["angle"])
    nat = {st: F(r1, f"R1_nativeA_{st}") for st in MESSAGES}
    band = {st: P1["blocks"][f"R1_nativeA_{st}"]["F"]["band"] for st in MESSAGES}
    inside = {st: band[st][0] <= nat[st] <= band[st][1] for st in MESSAGES}
    H.append({"id": "H1", "claim": "native F inside its band; fingerprint below the 0.05 rad floor (no N0 block)",
              "evidence": f"native A {nat['product_HT']:.3f} / {nat['bell_30_70']:.3f} (inside: {inside['product_HT']} / "
                          f"{inside['bell_30_70']}); fingerprint {fp:.3f} rad", "holds": all(inside.values()) and fp <= 0.05})
    rec = {st: (F(r2, f"R2_N1corr_{st}") - F(r1, f"R1_N1none_{st}")) / (nat[st] - F(r1, f"R1_N1none_{st}")) for st in MESSAGES}
    rs = {st: math.sqrt(S(r2, f"R2_N1corr_{st}") ** 2 + S(r1, f"R1_nativeA_{st}") ** 2) / (nat[st] - F(r1, f"R1_N1none_{st}"))
          for st in MESSAGES}
    H.append({"id": "H2", "claim": "TOMO-9 from real counts recovers the whole N1 loss (twin 1.00 +/- 0.06 / 1.00 +/- 0.02)",
              "evidence": "; ".join(f"{st}: {rec[st]:.3f} +/- {rs[st]:.3f}" for st in MESSAGES),
              "holds": all(abs(rec[st] - 1) <= max(2 * rs[st], 0.03) for st in MESSAGES), "values": rec})
    d = {st: F(r2, f"R2_N1corr_{st}") - F(r1, f"R1_N1oracle_{st}") for st in MESSAGES}
    ds = {st: math.sqrt(S(r2, f"R2_N1corr_{st}") ** 2 + S(r1, f"R1_N1oracle_{st}") ** 2) for st in MESSAGES}
    H.append({"id": "H3", "claim": "TOMO-9 = oracle within noise (TOMO-9 > oracle only with native coherent errors on the channel)",
              "evidence": "; ".join(f"{st}: TOMO-9 - oracle = {d[st]:+.4f} +/- {ds[st]:.4f}" for st in MESSAGES),
              "holds": all(abs(d[st]) <= 2 * ds[st] for st in MESSAGES), "values": d})
    n2 = {"none": F(r1, "R1_N2none_product_HT"), "tomo9": F(r2, "R2_N2corr_product_HT")}
    if "R2_DDN2none_product_HT" in r2["tomo"]:
        n2.update({"dd": F(r2, "R2_DDN2none_product_HT"), "tomo9_dd": F(r2, "R2_DDN2corr_product_HT")})
    H.append({"id": "H4", "claim": "N2: only the coherent part is recoverable; TOMO-9 alone ~ DD alone; TOMO-9 + DD over-corrects",
              "evidence": ", ".join(f"{k} {v:.3f}" for k, v in n2.items()), "values": n2, "holds": None})
    mit = {st: r1["tomo_mit"][f"R1_nativeA_{st}"]["F"] - nat[st] for st in MESSAGES}
    H.append({"id": "H5/H10", "claim": "readout biases the fidelity, not the correction (A); mitigated - raw ~ +0.01 in A",
              "evidence": "; ".join(f"{st}: mitigated - raw = {mit[st]:+.4f}" for st in MESSAGES), "values": mit,
              "holds": all(0.0 <= mit[st] <= 0.03 for st in MESSAGES)})
    mA, mB = r1["tomo_mit"]["R1_nativeA_product_HT"]["F"], r1["tomo_mit"]["R1_nativeB_product_HT"]["F"]
    mC = r1["tomo_mit"].get("R1_C_product_HT", {}).get("F")
    H.append({"id": "H6", "claim": "Bob-mitigated A > C > B, all within ~0.01 (product_HT)",
              "evidence": f"A {mA:.3f}, C {mC if mC is None else round(mC, 3)}, B {mB:.3f}",
              "holds": bool(mC is not None and mA >= mC >= mB and mA - mB <= 0.015)})
    keys = [k for k in r1["tomo"] if k in P1["blocks"]]
    below = sum(1 for k in keys if F(r1, k) < P1["blocks"][k]["F"]["prediction"])
    H.append({"id": "H7", "claim": "the twin is optimistic", "evidence": f"{below} of {len(keys)} Round 1 blocks below the L3 prediction",
              "holds": below > len(keys) / 2})
    pd_ = (r1["detector_fp"]["p"], r1["detector_cal"]["p"])
    H.append({"id": "H8", "claim": "the crosstalk detector is silent",
              "evidence": f"p = {pd_[0]:.3g} (fingerprint), {pd_[1]:.3g} (N1 calibration)", "holds": min(pd_) > 0.01})
    H.append({"id": "H9", "claim": "a native coherent rotation is claimed only above ~0.05 rad", "evidence": f"fingerprint {fp:.3f} rad",
              "holds": None, "native_rotation_claimed": fp > 0.05})
    gap = {st: r1["tomo_mit"][f"R1_nativeA_{st}"]["F"] - r1["tomo_mit"][f"R1_nativeB_{st}"]["F"] for st in MESSAGES}
    H.append({"id": "H11", "claim": "the A - B gap measures the feed-forward wait (twin: 0.011 / 0.011 mitigated)",
              "evidence": "; ".join(f"{st}: A - B (mitigated) = {gap[st]:+.4f}" for st in MESSAGES), "values": gap, "holds": None})
    return H


# --------------------------------------------------------------------------------------- recovered fraction
def fieller(N: float, Dn: float, vN: float, vD: float, cND: float, z: float) -> tuple[float, float]:
    """Fieller interval for f = N / D: the f with (N - f D)^2 <= z^2 (Var N - 2 f Cov(N, D) + f^2 Var D).
    Bounded when D^2 > z^2 Var D (true here)."""
    a = Dn ** 2 - z ** 2 * vD
    b = -2 * (N * Dn - z ** 2 * cND)
    c = N ** 2 - z ** 2 * vN
    disc = b * b - 4 * a * c
    lo, hi = sorted([(-b - math.sqrt(disc)) / (2 * a), (-b + math.sqrt(disc)) / (2 * a)])
    return lo, hi


FRACTION_BLOCKS = (("native", 1, "R1_nativeA"), ("none", 1, "R1_N1none"), ("corr", 2, "R2_N1corr"), ("oracle", 1, "R1_N1oracle"))


def recovered_fraction(data1: Mapping, data2: Mapping, n_boot: int = 2000, seed: int = 98501) -> dict:
    """f = (F_corr - F_none) / (F_native - F_none): the share of the injected N1 loss that the correction learned from
    real counts recovers, with the exact multinomial sigma of each linear F, the delta-method sigma, the 68 % and 95 %
    Fieller intervals and a multinomial bootstrap of the saved counts (``n_boot`` resamples, seeded)."""
    rows, boot = [], {}
    rng = np.random.default_rng(seed)
    for st in MESSAGES:
        w = linear_weights(message_state(st))
        c9 = {key: counts9(data1 if rnd == 1 else data2, prefix, st) for key, rnd, prefix in FRACTION_BLOCKS}
        Fv = {k: linear_F(c, w) for k, c in c9.items()}
        Sv = {k: math.sqrt(linear_F_variance(c, w)) for k, c in c9.items()}
        N, Dn = Fv["corr"] - Fv["none"], Fv["native"] - Fv["none"]
        f = N / Dn
        s_full = math.sqrt(Sv["corr"] ** 2 + f ** 2 * Sv["native"] ** 2 + (1 - f) ** 2 * Sv["none"] ** 2) / Dn
        vN, vD, cND = Sv["corr"] ** 2 + Sv["none"] ** 2, Sv["native"] ** 2 + Sv["none"] ** 2, Sv["none"] ** 2
        Fs = {}
        for key, _, _ in (("native", 1, ""), ("none", 1, ""), ("corr", 2, ""), ("oracle", 1, "")):
            acc = np.zeros(n_boot)
            for B in TOMO_BASES:
                n = np.array([c9[key][B].get(k, 0) for k in ("00", "01", "10", "11")], float)
                Nb = int(n.sum())
                acc += rng.multinomial(Nb, n / Nb, size=n_boot) @ w[B] / Nb
            Fs[key] = acc
        fb = (Fs["corr"] - Fs["none"]) / (Fs["native"] - Fs["none"])
        boot[st] = {"f_q": {str(q): float(np.quantile(fb, q)) for q in (0.025, 0.16, 0.5, 0.84, 0.975)},
                    "f_sd": float(np.std(fb, ddof=1)), "P_f_gt_1": float(np.mean(fb > 1)), "first10": fb[:10].tolist()}
        dor = Fv["corr"] - Fv["oracle"]
        rows.append({"state": st, "F": Fv, "sigma_exact": Sv, "f": f,
                     "sigma_8d": math.hypot(Sv["corr"], Sv["native"]) / Dn, "sigma_delta": s_full,
                     "fieller68": list(fieller(N, Dn, vN, vD, cND, 1.0)), "fieller95": list(fieller(N, Dn, vN, vD, cND, 1.96)),
                     "tomo_minus_oracle": dor, "tomo_minus_oracle_sd": math.hypot(Sv["corr"], Sv["oracle"])})
    return {"rows": rows, "bootstrap": {"seed": seed, "n_boot": n_boot, "results": boot}}


# --------------------------------------------------------------------------------------- job timing
def _iso(t):
    return dt.datetime.fromisoformat(str(t).replace("Z", "+00:00")) if t else None


def jobs_table() -> list[dict]:
    """One row per hardware job: shots, usage, pure circuit time, queue and run time (from the job records)."""
    rows = []
    for purpose, label in (("pilot", "pilot (R0)"), ("round1", "Round 1"), ("round2_main", "Round 2 main"),
                           ("round2_dd", "Round 2 DD")):
        raw = D.raw_job(purpose); rec = raw["record"]; m = raw.get("metrics") or {}
        ts = m.get("timestamps") or {}
        c, r, f = _iso(ts.get("created")), _iso(ts.get("running")), _iso(ts.get("finished"))
        est = [e["quantum_seconds"] for e in rec.get("ibm_estimates_queued", [])]
        row = {"job": label, "name": D.JOBS[purpose][1], "id": rec["id"], "shots": rec["shots_total"],
               "estimate_quick_s": rec["estimate_s"], "ibm_estimate_s": est[-1] if est else None, "usage_s": rec.get("usage_s"),
               "cap_s": rec["max_seconds"], "circuit_time_s": (m.get("circuits_execution_time_ns") or 0) * 1e-9,
               "queue_s": (r - c).total_seconds(), "run_wall_s": (f - r).total_seconds(),
               "created": ts.get("created"), "running": ts.get("running"), "finished": ts.get("finished")}
        row["ms_per_shot_measured"] = 1e3 * row["circuit_time_s"] / row["shots"]
        row["overhead_s"] = (row["usage_s"] or 0) - row["circuit_time_s"]
        rows.append(row)
    return rows


# --------------------------------------------------------------------------------------- everything
def reproduce(boot: int = 300, n_boot_fraction: int = 2000) -> dict:
    """The whole hardware analysis from the shipped raw counts (about 1 to 2 minutes)."""
    cfg = config()
    d1, d2 = D.load_round(1), D.load_round(2)
    r1 = analyse_round1(d1, cfg, boot=boot)
    corr = corrections_from_round1(r1, cfg)
    r2 = analyse_round2(d2, corr, r1["readout"]["M_bob"], cfg, boot=boot)
    P1, P2 = D.preregistration("predictions_round1.json"), D.preregistration("predictions_round2.json")
    A = D.preregistration("predictions_round2_addendum.json")
    v1, s1 = D.labels_info(1); v2, s2 = D.labels_info(2)
    chk1 = circuit_check(d1, P1["per_circuit"]["marginals"], P1["per_circuit"]["band"], s1, v1)
    chk2 = circuit_check(d2, A["per_circuit"]["marginals"], A["per_circuit"]["band"], s2, v2)
    tm = twirl_model(r1)
    cal = B.chain_calibration(D.calibration())
    eb = B.error_budget(r1, cal, tm["blocks"]["R1_nativeA_product_HT"]["model_F"], D.reference("8c_ladder.json"),
                        D.reference("8d_postdiction.json"), P1)
    return {"config": cfg, "r1": r1, "corr": corr, "r2": r2, "twirl_model": tm,
            "check_round1": chk1, "check_round2": chk2, "compare": compare(r1, r2, P1, P2, A),
            "hypotheses": hypotheses(r1, r2, P1), "recovered_fraction": recovered_fraction(d1, d2, n_boot_fraction),
            "chain_calibration": cal, "error_budget": eb, "n2": B.n2_analysis(r1, r2, cal),
            "dd": B.dd_scenarios(r1, r2, cal), "jobs": jobs_table()}
