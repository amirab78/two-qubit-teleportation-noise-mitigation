"""The replication of 25 October 2026 on ``ibm_kingston`` and ``ibm_fez``, from the raw counts.

Roadmap: step 8 (real hardware).

Eleven jobs, with the same estimators as the 4 October analysis (:mod:`.analysis`):

* ``K0`` / ``F0``: pilots (500 shots per circuit) against the twin's predicted marginals;
* ``K1`` / ``F1``: Round 1, three weeks later on Kingston and on a second device with no tuning: readout, the
  device alone, the injected N1 without correction, (Kingston) N1 + the correction stored from 4 October, the
  TOMO-9 of N1 and of the device fingerprint, (Kingston) waits of 4 / 8 / 16 us and a TOMO-9 at each wait;
* ``K2`` / ``F2``: Round 2, the fresh TOMO-9 corrections, a verification TOMO-9 and (Kingston) the corrected waits;
* ``K3`` / ``K4``: dynamical decoupling (DD) at 8 us with the idle message qubits pulsed / not pulsed: TOMO-9 under
  DD, DD alone, K1's correction + DD (Part 9.8's spectator test);
* ``G1`` (Fez): the device alone, N1 + F1's correction and the device + its fingerprint correction in one job;
* ``G2`` (Kingston): TOMO-9 under DD at 4 / 8 / 16 us (does DD's extra z rotation grow with the wait?);
* ``G3`` (Kingston): DD alone and DD + the correction learned under DD in G2, at 8 and 16 us.

Block keys: ``R<round>_<block>_<message>``, e.g. ``R2_N2corr16us_product_HT``. Round 5, 6, 7 are G1, G2, G3.
The pre-registered predictions are shipped (``preregistration/``); "model 2", the Pauli-twirl model from TOMO-9
matrices, is recomputed here from the counts; the twin ("model 1") is shipped as saved values only.
"""
from __future__ import annotations

import datetime as dt
import math
from typing import Mapping

import numpy as np

from .. import tomography as proc
from ..ceiling import procrustes
from ..detector import detect
from ..rotations import oracle_angles, so3, u_matrix
from ..states import MESSAGES, message_state
from ..twirl import leg_state_fidelity
from . import data as D
from .analysis import circuit_check, circuit_marginals, counts9, probe_dists
from .readout import READOUT_PATTERNS, confusion_from_counts
from .tomography import analyse as analyse_tomo

TAUS: tuple[int, ...] = (4000, 8000, 16000)                     # the waits (ns)
TAG = {4000: "4us", 8000: "8us", 16000: "16us"}
G3_TAUS: tuple[int, ...] = (8000, 16000)
DEVICE = {"K": "ibm_kingston", "F": "ibm_fez"}
#: Kingston's chain (logical q0..q5 -> physical qubit), as on 4 October
KINGSTON_LAYOUT: tuple[int, ...] = (41, 42, 36, 43, 21, 44)
#: plain names of the block keys (version A), without the round prefix and the wait
NAMES = {"nativeA": "device alone", "N1none": "N1, no correction", "N1old": "N1, correction stored from 4 Oct",
         "N1corr": "N1, TOMO-9 correction", "native": "device alone", "fpcorr": "device alone, fingerprint correction",
         "N2none": "wait, no correction", "N2corr": "wait, TOMO-9 correction", "DDnone": "DD alone, message pulsed",
         "DDcorr": "K1 correction + DD, message pulsed", "DDUnone": "DD alone",
         "DDUcorr": "correction + DD"}


def describe(key: str) -> tuple[str, str]:
    """``"R2_N2corr16us_product_HT"`` -> ("16 µs wait, TOMO-9 correction", "product_HT"): plain name and message."""
    st = "bell_30_70" if key.endswith("bell_30_70") else "product_HT"
    block = key[: -len(st) - 1].split("_", 1)[1]
    wait = next((t for t in ("16us", "4us", "8us") if t in block), None)
    base = block.replace(wait, "") if wait else block
    if base.startswith("DD") and not wait:
        wait = "8us"
    name = NAMES.get(base, base)
    if base == "DDUcorr":
        name = "K1 correction + DD" if key.startswith("R4_") else "correction learned under DD + DD"
    if wait:
        w = wait.replace("us", " µs")
        name = f"{w} wait, " + name.replace("wait, ", "") if base.startswith("N2") else f"{name} ({w})"
    return name, st


# --------------------------------------------------------------------------------------- the estimators
def analyse(data: Mapping, shots: Mapping[str, int], boot: int = 300, M=None) -> dict:
    """Every tomography block (raw and Bob-mitigated), the readout calibration when the job has one (its Bob matrix
    is then used), every TOMO-9 estimate and the crosstalk detector, as in the 4 October analysis."""
    res: dict = {"tomo": {}, "tomo_mit": {}, "tomo9": {}, "detector": {}}
    if all(f"readout_{p}" in data for p in READOUT_PATTERNS):
        res["readout"] = confusion_from_counts({p: data[f"readout_{p}"]["counts"] for p in READOUT_PATTERNS})
        M = res["readout"]["M_bob"]
    res["M_used"] = M
    blocks, t9 = {}, {}
    for lab in data:
        if "_tomo9_" in lab:
            t9.setdefault(lab.split("_tomo9_")[0], shots[lab])
        elif lab[:1] == "R" and lab[1:2].isdigit():
            for st in MESSAGES:
                if f"_{st}_" in lab:
                    blocks[(lab.split(f"_{st}_")[0], st)] = True
    for pre, st in blocks:
        psi = message_state(st)
        res["tomo"][f"{pre}_{st}"] = _strip(analyse_tomo(counts9(data, pre, st, "A"), psi, boot))
        if M is not None:
            res["tomo_mit"][f"{pre}_{st}"] = _strip(analyse_tomo(counts9(data, pre, st, "A", M), psi, 0))
    lab9 = proc.probe_labels("tomo9")
    for pre, n in t9.items():
        d = probe_dists(data, [f"{pre}_{l}" for l in lab9])
        key = pre.split("_", 1)[1]
        res["tomo9"][key] = proc.estimate(d, "tomo9")
        res["detector"][key] = detect(d, n)
    return res


def _strip(a: dict) -> dict:
    return {k: v for k, v in a.items() if k not in ("rho_real", "rho_imag")}


def _so3_diff(x, y) -> float:
    """Largest entry difference of the two legs' Bloch rotations of two corrections (6 U angles each)."""
    return max(float(np.abs(so3(u_matrix(*x[3 * i:3 * i + 3])) - so3(u_matrix(*y[3 * i:3 * i + 3]))).max())
               for i in range(2))


def _zrot(est: Mapping) -> dict[int, float]:
    """z component of the rotation vector TOMO-9 finds on each leg."""
    return {q: est[f"rotation_{q}"]["angle"] * est[f"rotation_{q}"]["axis"][2] for q in (2, 3)}


def zboot(data: Mapping, prefix: str, boot: int = 300, seed: int = 2560) -> dict:
    """z part of the TOMO-9 rotation of one calibration (label prefix), with a multinomial bootstrap sigma."""
    dists = probe_dists(data, [f"{prefix}_{l}" for l in proc.probe_labels("tomo9")])
    est = proc.estimate(dists, "tomo9")
    rng = np.random.default_rng(seed)
    zs = []
    for _ in range(boot):
        res = []
        for d in dists:
            ks = list(d); n = int(sum(d.values())); p = np.array([d[k] for k in ks], float) / n
            res.append(dict(zip(ks, rng.multinomial(n, p))))
        z = _zrot(proc.estimate(res, "tomo9")); zs.append([z[2], z[3]])
    zs = np.array(zs)
    return {"z": _zrot(est), "sigma": {q: float(zs[:, i].std()) for i, q in enumerate((2, 3))}, "est": est}


def _check(job: str, data: Mapping, shots: Mapping, versions: Mapping, pc: Mapping) -> list[dict]:
    return circuit_check(data, pc["marginals"], pc["band"], shots, versions)


def _prediction(job: str) -> dict:
    return D.replication_preregistration(D.replication_frozen(job)["prediction_file"])


# --------------------------------------------------------------------------------------- model 2 (twirl model)
def so3_of_x(x6) -> dict[int, np.ndarray]:
    """The two legs' Bloch rotations of a correction given as 6 U angles."""
    return {2: so3(u_matrix(*x6[:3])), 3: so3(u_matrix(*x6[3:]))}


def twirl_F(st: str, M2, M3) -> float:
    """Bob's fidelity in the twirl model: only diag(M) of each leg reaches Bob."""
    return leg_state_fidelity(message_state(st), np.diag(np.diag(M2)), np.diag(np.diag(M3)))


def T_of(est: Mapping) -> dict[int, np.ndarray]:
    t = est["T"]
    return {q: np.array(t[str(q)] if str(q) in t else t[q], float) for q in (2, 3)}


def _rz3(a: float) -> np.ndarray:
    c, s = math.cos(a), math.sin(a)
    return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1.0]])


def spectator_rad(props: Mapping, tau_us: float = 8.0) -> dict[int, float]:
    """pi zeta t of the spectator couplers (message qubit q0 next to q2, q1 next to q3) from a calibration record."""
    zz = {x["name"]: x["value"] for x in props["general"] if x["name"].startswith("zz_")}
    lay = KINGSTON_LAYOUT
    z = lambda a, b: zz.get(f"zz_{a}{b}", zz.get(f"zz_{b}{a}")) * 1e6           # noqa: E731  GHz -> kHz
    return {2: math.pi * z(lay[2], lay[0]) * 1e3 * tau_us * 1e-6, 3: math.pi * z(lay[3], lay[1]) * 1e3 * tau_us * 1e-6}


def dd_model(r1: Mapping, corr8, props: Mapping) -> dict:
    """Part 9.8's DD model at 8 us from K1's TOMO-9 matrices: T(N2) = Rc . D; DD refocuses every z phase and keeps the
    non-z native rotation Rf (from the fingerprint); pulsed spectators add R_z(+pi zeta t)."""
    from scipy.spatial.transform import Rotation
    Tn8, Tfp = T_of(r1["tomo9"]["calN28us"]), T_of(r1["tomo9"]["fp"])
    Rc = {q: procrustes(Tn8[q]).T for q in (2, 3)}
    D8 = {q: Rc[q].T @ Tn8[q] for q in (2, 3)}
    Rf = {q: procrustes(Tfp[q]).T for q in (2, 3)}
    Rcorr = so3_of_x(corr8)
    s = spectator_rad(props, 8.0)
    out = {}
    for opt, pulsed, tg, rn in (("dd", True, "DD", 3), ("ddu", False, "DDU", 4)):
        kept = {}
        for q in (2, 3):
            vv = Rotation.from_matrix(Rf[q]).as_rotvec(); vv[2] = 0.0
            nonz = Rotation.from_rotvec(vv).as_matrix()
            kept[q] = (_rz3(s[q]) @ nonz) if pulsed else nonz
        Tdd = {q: kept[q] @ D8[q] for q in (2, 3)}
        rv = {q: Rotation.from_matrix(kept[q]).as_rotvec().tolist() for q in (2, 3)}
        out[opt] = {"pulsed_spectators": pulsed, "spectator_rad": s, "kept_rotvec": rv,
                    "kept_z_rad": {q: rv[q][2] for q in (2, 3)},
                    f"R{rn}_{tg}none_product_HT": twirl_F("product_HT", Tdd[2], Tdd[3]),
                    f"R{rn}_{tg}corr_product_HT": twirl_F("product_HT", Rcorr[2] @ Tdd[2], Rcorr[3] @ Tdd[3])}
    return out


def model2_round2(r1: Mapping, corr: Mapping, kingston: bool) -> dict:
    """Model 2 of the Round 2 blocks from Round 1's own TOMO-9 matrices and the learned corrections."""
    Tc9 = T_of(r1["tomo9"]["cal9"]); Rn1 = so3_of_x(corr["N1"])
    out = {f"R2_N1corr_{st}": twirl_F(st, Rn1[2] @ Tc9[2], Rn1[3] @ Tc9[3]) for st in MESSAGES}
    for st in MESSAGES:
        out[f"R1_nativeA_{st} (measured)"] = r1["tomo"][f"R1_nativeA_{st}"]["F"]
    for tau in TAUS if kingston else ():
        t = TAG[tau]; Tn = T_of(r1["tomo9"][f"calN2{t}"]); Rc = so3_of_x(corr[f"N2_{t}"])
        out[f"R2_N2corr{t}_product_HT"] = twirl_F("product_HT", Rc[2] @ Tn[2], Rc[3] @ Tn[3])
    return out


def wait_khz(r1: Mapping) -> dict:
    """The waiting detuning from K1: least-squares slope through the origin of the TOMO-9 rotation angle at 4, 8 and
    16 us, as a frequency (angle / (2 pi t)), with each wait's angle and the fit residuals."""
    t = np.array([tau / 1000.0 for tau in TAUS])
    out = {"t_us": t.tolist()}
    for q in (2, 3):
        th = np.array([r1["tomo9"][f"calN2{TAG[tau]}"][f"rotation_{q}"]["angle"] for tau in TAUS])
        slope = float(th @ t / (t @ t))
        out[f"q{q}"] = slope / (2 * np.pi) * 1e3
        out[f"angles_q{q}"] = th.tolist()
        out[f"residuals_q{q}"] = (th - slope * t).tolist()
    return out


# --------------------------------------------------------------------------------------- the jobs
def pilot(job: str) -> dict:
    """K0 / F0: the per-circuit check against the twin and the physics of every pilot circuit."""
    data, shots, vers = D.replication_job(job)
    P = _prediction(job)["pilot"]
    chk = circuit_check(data, P["marginals"], P["band"], shots, vers)
    phys = []
    for lab, e in data.items():
        m = circuit_marginals(e, lab, vers[lab]); p = P["marginals"][lab]
        if "_readout_" in lab:
            i = 0 if lab.endswith("000000") else 15
            what, a, b = "P(all 4 Alice bits = pattern)", m["alice"][i], p["alice"][i]
        else:
            tail = lab.split("_")[-1]
            what = f"Bob <{tail}> parity (frame-corrected)"
            a, b = _parity(m["bob"]), _parity(p["bob"])
        phys.append({"label": lab, "quantity": what, "hardware": float(a), "twin": float(b), "diff": float(a - b)})
    return {"check": chk, "n_flagged": sum(1 for r in chk if not r["ok"]), "physics": phys}


def _parity(p) -> float:
    p = np.asarray(p, float)
    return float(p[0] - p[1] - p[2] + p[3])


def corrections(r1: Mapping) -> dict:
    out = {"N1": list(r1["tomo9"]["cal9"]["x_new"])}
    for tau in TAUS:
        k = f"calN2{TAG[tau]}"
        if k in r1["tomo9"]:
            out[f"N2_{TAG[tau]}"] = list(r1["tomo9"][k]["x_new"])
    return out


def round1(job: str, boot: int = 300) -> dict:
    """K1 / F1: every block, the corrections for Round 2, the per-circuit check and the pre-registered go / no-go rules."""
    data, shots, vers = D.replication_job(job)
    r1 = analyse(data, shots, boot)
    P = _prediction(job)
    chk = _check(job, data, shots, vers, P["per_circuit"])
    corr = corrections(r1)
    lim = P["limits"]
    rules = []

    def rule(n, ok, v, l, act):
        rules.append({"rule": n, "ok": bool(ok), "value": v, "limit": l, "if_failed": act})
    nf = sum(1 for r in chk if not r["ok"])
    rule("per-circuit check: no circuit far from its prediction", nf == 0, nf, 0, "stop and inspect the flagged circuits")
    emax = max(r1["readout"]["e0"] + r1["readout"]["e1"])
    rule("readout errors plausible", emax < lim["readout_max"], emax, lim["readout_max"], "stop; the chain's readout is broken today")
    for st in MESSAGES:
        F = r1["tomo"][f"R1_nativeA_{st}"]["F"]
        rule(f"native fidelity ({st}) not far below the prediction", F >= lim[f"native_A_min_{st}"], F, lim[f"native_A_min_{st}"],
             "stop; the device is much worse than its calibration")
    dist = _so3_diff(corr["N1"], oracle_angles("N1"))
    rule("TOMO-9 correction of N1 close to the exact inverse", dist <= lim["cal9_oracle_max"], dist, lim["cal9_oracle_max"],
         "stop; the calibration is wrong")
    gain = r1["tomo9"]["cal9"]["Favg_after"] - r1["tomo9"]["cal9"]["Favg_before"]
    rule("TOMO-9 predicts a large gain on N1", gain >= lim["cal9_gain_min"], gain, lim["cal9_gain_min"], "stop; N1 did not act as designed")
    return {"r1": r1, "corr": corr, "check": chk, "n_flagged": nf, "go_no_go": rules}


def round2(job: str, R1: Mapping, boot: int = 300) -> dict:
    """K2 / F2: the fresh corrections (Bob mitigated with Round 1's readout), model 2 recomputed from Round 1's counts,
    the per-circuit check and the pre-registered REPLICATION and WAIT criteria."""
    kingston = job[0] == "K"
    data, shots, vers = D.replication_job(job)
    r2 = analyse(data, shots, boot, M=R1["r1"]["readout"]["M_bob"])
    P = _prediction(job)
    chk = _check(job, data, shots, vers, P["per_circuit"]["K2"])
    m2 = model2_round2(R1["r1"], R1["corr"], kingston)
    t1, t2 = R1["r1"]["tomo"], r2["tomo"]
    blocks, crit = [], []
    for st in MESSAGES:
        k = f"R2_N1corr_{st}"; F, s = t2[k]["F"], t2[k]["F_sigma"]
        Fn, sn = t1[f"R1_nativeA_{st}"]["F"], t1[f"R1_nativeA_{st}"]["F_sigma"]
        F0 = t1[f"R1_N1none_{st}"]["F"]
        frac = (F - F0) / (Fn - F0)
        blocks.append({"block": k, "F": F, "sigma": s, "F_mit": r2["tomo_mit"][k]["F"], "model_2": m2[k],
                       "twin": P["twin_K2"]["blocks"][k]["prediction"], "z_vs_model_2": (F - m2[k]) / s,
                       "reference": f"native {Fn:.3f}, none {F0:.3f}", "recovered_fraction": frac,
                       "stored_4oct_correction": t1.get(f"R1_N1old_{st}", {}).get("F")})
        crit.append({"criterion": f"REPLICATION ({st}): F >= F(native) - 2 sigma", "ok": F >= Fn - 2 * np.hypot(s, sn),
                     "value": F, "limit": Fn - 2 * np.hypot(s, sn)})
        crit.append({"criterion": f"REPLICATION ({st}): recovered fraction >= 0.9", "ok": frac >= 0.9, "value": frac, "limit": 0.9})
    for tau in TAUS if kingston else ():
        t = TAG[tau]; k = f"R2_N2corr{t}_product_HT"; F, s = t2[k]["F"], t2[k]["F_sigma"]
        F0, s0 = t1[f"R1_N2none{t}_product_HT"]["F"], t1[f"R1_N2none{t}_product_HT"]["F_sigma"]
        blocks.append({"block": k, "F": F, "sigma": s, "F_mit": r2["tomo_mit"][k]["F"], "model_2": m2[k],
                       "twin": P["twin_K2"]["blocks"][k]["prediction"], "z_vs_model_2": (F - m2[k]) / s,
                       "reference": f"none {F0:.3f}", "gain": F - F0, "gain_sigma": float(np.hypot(s, s0))})
        crit.append({"criterion": f"WAIT {t}: F(corrected) - F(none) > 3 sigma", "ok": F - F0 > 3 * np.hypot(s, s0),
                     "value": (F - F0) / np.hypot(s, s0), "limit": 3.0})
    return {"r2": r2, "blocks": blocks, "criteria": crit, "check": chk, "n_flagged": sum(1 for r in chk if not r["ok"]),
            "verification": r2["tomo9"].get("ver"), "detector": r2["detector"], "model_2": m2}


def dd(job: str, K1: Mapping, K2: Mapping, boot: int = 300, seed: int = 2540) -> dict:
    """K3 (message qubits pulsed) or K4 (not pulsed): DD alone and K1's correction + DD, the rotation TOMO-9 finds
    under DD (z part, bootstrap sigma) against Part 9.8's model (recomputed from K1) and the twin."""
    rnd, tag = (3, "DD") if job == "K3" else (4, "DDU")
    opt = "dd" if job == "K3" else "ddu"
    data, shots, vers = D.replication_job(job)
    r = analyse(data, shots, boot, M=K1["r1"]["readout"]["M_bob"])
    P = _prediction(job)
    chk = _check(job, data, shots, vers, P["per_circuit"][job])
    key = f"{tag}cal"
    Z = zboot(data, f"R{rnd}_{key}", boot, seed)
    est = r["tomo9"][key]
    props = D.replication_calibration("ibm_kingston", K1_STAMP)
    m2 = dd_model(K1["r1"], K1["corr"]["N2_8us"], props)[opt]
    rot = [{"leg": f"q{q}", "z measured (rad)": Z["z"][q], "sigma": Z["sigma"][q], "z model 2 (rad)": m2["kept_z_rad"][q],
            "pull": (Z["z"][q] - m2["kept_z_rad"][q]) / Z["sigma"][q], "angle": est[f"rotation_{q}"]["angle"],
            "axis": est[f"rotation_{q}"]["axis"], "shrink": est[f"shrink_{q}"]} for q in (2, 3)]
    blocks = []
    for k in (f"R{rnd}_{tag}none_product_HT", f"R{rnd}_{tag}corr_product_HT"):
        t = r["tomo"][k]
        blocks.append({"block": k, "F": t["F"], "sigma": t["F_sigma"], "F_mit": r["tomo_mit"][k]["F"], "model_2": m2[k],
                       "twin": P["twin_DD"][job]["blocks"][k]["prediction"], "z_vs_model_2": (t["F"] - m2[k]) / t["F_sigma"]})
    t1 = K1["r1"]["tomo"]
    return {"pulsed_spectators": opt == "dd", "blocks": blocks, "rotation": rot, "tomo9": est, "detector": r["detector"],
            "check": chk, "n_flagged": sum(1 for c in chk if not c["ok"]), "model_2": m2,
            "reference_K1": {"N2none8us": t1["R1_N2none8us_product_HT"]["F"], "native": t1["R1_nativeA_product_HT"]["F"]},
            "reference_K2": {"N2corr8us": next(b["F"] for b in K2["blocks"] if b["block"] == "R2_N2corr8us_product_HT")}}


K1_STAMP = "20261025T0641Z"            # the calibration record in force for K0 and K1


def g1(F1: Mapping, F2: Mapping, boot: int = 300) -> dict:
    """G1 (Fez): the device alone, N1 + F1's correction and the device + its fingerprint correction in ONE job: is F2's
    'corrected above native' real (REAL = the two-job gap F2 - F1) or drift between jobs (DRIFT = model 2's gap)?"""
    data, shots, vers = D.replication_job("G1")
    r = analyse(data, shots, boot, M=F1["r1"]["readout"]["M_bob"])
    P = _prediction("G1")["G1"]
    chk = _check("G1", data, shots, vers, P["per_circuit"])
    c1, cfp = F1["corr"]["N1"], list(F1["r1"]["tomo9"]["fp"]["x_new"])
    Tfp, Tc9 = T_of(F1["r1"]["tomo9"]["fp"]), T_of(F1["r1"]["tomo9"]["cal9"])
    Rn1, Rfp = so3_of_x(c1), so3_of_x(cfp)
    m2 = {}
    for st in MESSAGES:
        m2[f"R5_native_{st}"] = twirl_F(st, Tfp[2], Tfp[3])
        m2[f"R5_N1corr_{st}"] = twirl_F(st, Rn1[2] @ Tc9[2], Rn1[3] @ Tc9[3])
        m2[f"R5_fpcorr_{st}"] = twirl_F(st, Rfp[2] @ Tfp[2], Rfp[3] @ Tfp[3])
    f2 = {b["block"].replace("R2_N1corr_", ""): b["F"] for b in F2["blocks"]}
    hyp = {st: {"REAL (F2 - F1, two jobs)": f2[st] - F1["r1"]["tomo"][f"R1_nativeA_{st}"]["F"],
                "DRIFT (model 2)": m2[f"R5_N1corr_{st}"] - m2[f"R5_native_{st}"]} for st in MESSAGES}
    t = r["tomo"]
    blocks = [{"block": k, "F": v["F"], "sigma": v["F_sigma"], "F_mit": r["tomo_mit"][k]["F"], "model_2": m2[k],
               "twin": P["twin"][k]["prediction"]} for k, v in t.items()]
    tests = []
    for st in MESSAGES:
        for what, k in (("N1 + correction", f"R5_N1corr_{st}"), ("fingerprint correction", f"R5_fpcorr_{st}")):
            Dv = t[k]["F"] - t[f"R5_native_{st}"]["F"]; s = float(np.hypot(t[k]["F_sigma"], t[f"R5_native_{st}"]["F_sigma"]))
            e = {"state": st, "comparison": f"{what} - native (same job)", "D": Dv, "sigma": s, "z": Dv / s}
            if what == "N1 + correction":
                real, drift = hyp[st]["REAL (F2 - F1, two jobs)"], hyp[st]["DRIFT (model 2)"]
                e.update({"REAL": real, "DRIFT": drift,
                          "log10 LR (REAL / DRIFT)": float(((Dv - drift) ** 2 - (Dv - real) ** 2) / (2 * s * s) / np.log(10))})
            tests.append(e)
    return {"blocks": blocks, "tests": tests, "check": chk, "n_flagged": sum(1 for c in chk if not c["ok"]),
            "model_2": m2, "gap_hypotheses": hyp}


def g2(boot: int = 300) -> dict:
    """G2 (Kingston): TOMO-9 under DD (message qubits not pulsed) at 4, 8 and 16 us; z(t) = a + b t on each leg."""
    data, shots, vers = D.replication_job("G2")
    P = _prediction("G2")["G2"]
    chk = _check("G2", data, shots, vers, P["per_circuit"])
    t = np.array([tau / 1000.0 for tau in TAUS])
    Z = {TAG[tau]: zboot(data, f"R6_DDU{TAG[tau]}cal", boot, 2563 + i) for i, tau in enumerate(TAUS)}
    rot, fits = [], {}
    for q in (2, 3):
        z = np.array([Z[TAG[tau]]["z"][q] for tau in TAUS]); s = np.array([Z[TAG[tau]]["sigma"][q] for tau in TAUS])
        A = np.vstack([np.ones_like(t), t]).T / s[:, None]
        cov = np.linalg.inv(A.T @ A); ab = cov @ A.T @ (z / s)
        fits[f"q{q}"] = {"a_rad": float(ab[0]), "a_sigma": float(np.sqrt(cov[0, 0])), "b_rad_per_us": float(ab[1]),
                         "b_sigma": float(np.sqrt(cov[1, 1])), "chi2_1dof": float(np.sum(((z - (ab[0] + ab[1] * t)) / s) ** 2)),
                         "CONST favoured (b ~ 0)": bool(abs(ab[1]) < 2 * np.sqrt(cov[1, 1])),
                         "LINEAR favoured (a ~ 0)": bool(abs(ab[0]) < 2 * np.sqrt(cov[0, 0]))}
        for tau in TAUS:
            tg = TAG[tau]; h = P["z_hypotheses"][tg]
            rot.append({"wait": tg, "leg": f"q{q}", "z measured (rad)": Z[tg]["z"][q], "sigma": Z[tg]["sigma"][q],
                        "CONST": h["CONST"][str(q)], "LINEAR": h["LINEAR"][str(q)],
                        "shrink": Z[tg]["est"][f"shrink_{q}"], "angle": Z[tg]["est"][f"rotation_{q}"]["angle"]})
    corr = {TAG[tau]: list(Z[TAG[tau]]["est"]["x_new"]) for tau in G3_TAUS}
    est = {TAG[tau]: Z[TAG[tau]]["est"] for tau in TAUS}
    return {"rotation": rot, "fits": fits, "check": chk, "n_flagged": sum(1 for c in chk if not c["ok"]),
            "corrections_for_G3": corr, "tomo9": est}


def g3(G2: Mapping, K1: Mapping, K2: Mapping, K4: Mapping, boot: int = 300) -> dict:
    """G3 (Kingston): DD alone and DD + G2's correction (learned under DD) at 8 and 16 us; model 2 from G2's TOMO-9
    matrices; G3-A (the correction helps) and G3-B (against K2's fresh correction without DD)."""
    data, shots, vers = D.replication_job("G3")
    P = _prediction("G3")
    corr = G2["corrections_for_G3"]
    r = analyse(data, shots, boot, M=K1["r1"]["readout"]["M_bob"])
    chk = _check("G3", data, shots, vers, P["per_circuit"])
    m2 = {}
    for tau in G3_TAUS:
        t = TAG[tau]; T = T_of(G2["tomo9"][t]); Rc = so3_of_x(corr[t])
        m2[f"R7_DDU{t}none_product_HT"] = twirl_F("product_HT", T[2], T[3])
        m2[f"R7_DDU{t}corr_product_HT"] = twirl_F("product_HT", Rc[2] @ T[2], Rc[3] @ T[3])
    k2 = {b["block"]: b for b in K2["blocks"]}
    ref = {t: {"K2 fresh correction, no DD": [k2[f"R2_N2corr{t}_product_HT"]["F"], k2[f"R2_N2corr{t}_product_HT"]["sigma"]]}
           for t in ("8us", "16us")}
    k4 = {b["block"]: b for b in K4["blocks"]}
    ref["8us"]["K4 DD alone"] = [k4["R4_DDUnone_product_HT"]["F"], k4["R4_DDUnone_product_HT"]["sigma"]]
    blocks = [{"block": k, "F": v["F"], "sigma": v["F_sigma"], "F_mit": r["tomo_mit"][k]["F"], "model_2": m2[k],
               "twin": P["twin"][k]["prediction"], "pull_model_2": (v["F"] - m2[k]) / v["F_sigma"]} for k, v in r["tomo"].items()]
    crit = []
    for tau in G3_TAUS:
        t = TAG[tau]; a, b = r["tomo"][f"R7_DDU{t}corr_product_HT"], r["tomo"][f"R7_DDU{t}none_product_HT"]
        s = float(np.hypot(a["F_sigma"], b["F_sigma"]))
        crit.append({"criterion": f"G3-A {t}: F(DD + corr) - F(DD alone) > 3 sigma", "value": a["F"] - b["F"], "sigma": s,
                     "z": (a["F"] - b["F"]) / s, "ok": a["F"] - b["F"] > 3 * s})
        kf, ks = ref[t]["K2 fresh correction, no DD"]
        s2 = float(np.hypot(a["F_sigma"], ks)); dz = (a["F"] - kf) / s2
        crit.append({"criterion": f"G3-B {t}: DD + corr vs K2 fresh correction without DD ({kf:.3f})", "value": a["F"] - kf,
                     "sigma": s2, "z": dz, "ok": None, "reading": "DD + correction better" if dz > 2 else ("worse" if dz < -2 else "equal")})
    return {"blocks": blocks, "criteria": crit, "check": chk, "n_flagged": sum(1 for c in chk if not c["ok"]),
            "model_2": m2, "references": ref}


# --------------------------------------------------------------------------------------- job timing
def _iso(t):
    return dt.datetime.fromisoformat(str(t).replace("Z", "+00:00")) if t else None


def jobs_table() -> list[dict]:
    """One row per replication job: device, shots, IBM's detailed estimate, cap, usage, circuit time, queue, run time."""
    rows = []
    for job, (jid, dev) in D.REPLICATION_JOBS.items():
        raw = D.replication_raw(job); rec = raw["record"]; man = D.replication_frozen(job)
        ts = (raw.get("metrics") or {}).get("timestamps") or {}
        c, r, f = _iso(ts.get("created")), _iso(ts.get("running")), _iso(ts.get("finished"))
        rows.append({"job": job, "device": dev, "id": jid, "circuits": man["n_bindings"], "shots": rec["shots_total"],
                     "dd": rec["dd"], "dd_options": rec.get("x_options"), "estimate_detailed_s": man["estimate_detailed_s"],
                     "ibm_estimate_s": rec.get("ibm_estimate_s"), "cap_s": rec["max_seconds"], "usage_s": rec["usage_s"],
                     "circuit_time_s": (raw["metrics"].get("circuits_execution_time_ns") or 0) * 1e-9,
                     "queue_s": (r - c).total_seconds(), "run_wall_s": (f - r).total_seconds(), "created": ts.get("created")})
    return rows


# --------------------------------------------------------------------------------------- everything
def reproduce(boot: int = 300) -> dict:
    """The whole replication analysis from the shipped raw counts (a few minutes)."""
    out: dict = {"K0": pilot("K0"), "F0": pilot("F0")}
    out["K1"] = round1("K1", boot); out["F1"] = round1("F1", boot)
    out["K2"] = round2("K2", out["K1"], boot); out["F2"] = round2("F2", out["F1"], boot)
    out["K3"] = dd("K3", out["K1"], out["K2"], boot); out["K4"] = dd("K4", out["K1"], out["K2"], boot)
    out["G1"] = g1(out["F1"], out["F2"], boot)
    out["G2"] = g2(boot)
    out["G3"] = g3(out["G2"], out["K1"], out["K2"], out["K4"], boot)
    out["wait_khz"] = wait_khz(out["K1"]["r1"])
    out["jobs"] = jobs_table()
    return out
