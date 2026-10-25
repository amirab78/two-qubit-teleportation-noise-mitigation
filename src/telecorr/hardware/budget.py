"""Where the hardware fidelity goes: the readout share of the native error, and the 8 us waiting channel (N2).

Roadmap: step 8 (real hardware).

Both analyses use only the Round 1 TOMO-9 results, the readout measured in Round 1 and the device calibration
record in force when the jobs ran (T1, T2, readout and static ZZ of the six qubits of the chain).
"""
from __future__ import annotations

import math
from typing import Mapping

import numpy as np

from ..ceiling import procrustes
from ..states import message_state
from ..twirl import leg_state_fidelity

#: physical qubits of the chain on ibm_kingston for the logical roles q0..q5 (chosen in Stage 8a)
CHAIN: dict[str, int] = {"q0": 41, "q1": 42, "q2": 36, "q3": 43, "q4": 21, "q5": 44}
CHAIN_ORDER: tuple[int, ...] = (21, 36, 41, 42, 43, 44)                 # path order q4 - q2 - q0 - q1 - q3 - q5
TAU_US: float = 8.0                                                    # the N2 wait


def chain_calibration(props: Mapping) -> dict:
    """T1, T2 (us) and readout errors of the chain's qubits, and the static ZZ (kHz) of its couplers, from a
    calibration record (BackendProperties as JSON)."""
    qubits = {}
    for role, p in CHAIN.items():
        v = {x["name"]: x["value"] for x in props["qubits"][p]}
        qubits[role] = {"phys": p, **{k: v.get(k) for k in ("T1", "T2", "readout_error", "prob_meas1_prep0", "prob_meas0_prep1")}}
    zz = {x["name"]: x["value"] for x in props["general"] if x["name"].startswith("zz_")}
    edges = {}
    for a, b in zip(CHAIN_ORDER[:-1], CHAIN_ORDER[1:]):
        val = zz.get(f"zz_{a}{b}", zz.get(f"zz_{b}{a}"))
        edges[f"{a}-{b}"] = {"zz_kHz": None if val is None else val * 1e6}      # GHz -> kHz
    return {"qubits": qubits, "edges": edges}


def leg_F_product(p) -> float:
    """Fidelity of one product_HT leg (Bloch vector in the xy plane at 45 degrees) through the Pauli channel
    p = [pI, pX, pY, pZ]: F_leg = 1 - pZ - (pX + pY) / 2."""
    return 1 - p[3] - (p[1] + p[2]) / 2


def _readout_pauli(e) -> dict[int, list[float]]:
    """Bob's Pauli errors caused by readout alone: Alice's message bit (q0 / q1) decides Bob's Z correction, her
    channel bit (q2 / q3) his X correction, and Bob's own readout acts on every basis (about e / 2 per Pauli)."""
    P = {}
    for q, (a, c, b) in ((2, (0, 2, 4)), (3, (1, 3, 5))):
        P[q] = [0.0, e[c] + e[b] / 2, e[b] / 2, e[a] + e[b] / 2]
        P[q][0] = 1 - sum(P[q][1:])
    return P


def error_budget(r1: Mapping, cal: Mapping, model_native_product: float, ladder: Mapping | None = None,
                 postdiction: Mapping | None = None, P1: Mapping | None = None) -> dict:
    """Readout share of the native error (measured readout vs the calibration record), the Pauli fingerprint of each
    learned channel, and (with the saved twin results) the waterfall ideal -> twin -> hardware."""
    ro = r1["readout"]
    eb = [(ro["e0"][i] + ro["e1"][i]) / 2 for i in range(6)]
    rec = cal["qubits"]
    eb_rec = [(rec[f"q{i}"]["prob_meas1_prep0"] + rec[f"q{i}"]["prob_meas0_prep1"]) / 2 for i in range(6)]
    legs: dict = {}
    for kind in ("fingerprint", "cal9", "calN2"):
        legs[kind] = {q: {"pauli": r1[kind][f"pauli_{q}"], "F_leg_product": leg_F_product(r1[kind][f"pauli_{q}"])} for q in (2, 3)}
        legs[kind]["F_product_formula"] = legs[kind][2]["F_leg_product"] * legs[kind][3]["F_leg_product"]
    ro_P, rec_P = _readout_pauli(eb), _readout_pauli(eb_rec)
    out = {"readout_bit_flip": {"measured": eb, "record": eb_rec}, "legs": legs,
           "readout_share": {2: {"Z (Alice q0 = Q41)": eb[0], "X (Alice q2 = Q36)": eb[2], "X, Y, Z each (Bob q4 = Q21)": eb[4] / 2},
                             3: {"Z (Alice q1 = Q42)": eb[1], "X (Alice q3 = Q43)": eb[3], "X, Y, Z each (Bob q5 = Q44)": eb[5] / 2}},
           "readout_pauli": ro_P, "readout_pauli_record": rec_P,
           "F_product_readout_only": leg_F_product(ro_P[2]) * leg_F_product(ro_P[3]),
           "F_product_readout_only_record": leg_F_product(rec_P[2]) * leg_F_product(rec_P[3]),
           "model_native_product": model_native_product,
           "hardware_native_product": r1["tomo"]["R1_nativeA_product_HT"]["F"],
           "q41_cost_product": eb[0] * legs["fingerprint"][3]["F_leg_product"]}
    if ladder and postdiction and P1:
        lad = ladder["levels"]
        water = {}
        for st, la in (("product_HT", "A_product_HT"), ("bell_30_70", "A_bell_30_70")):
            k = f"R1_nativeA_{st}"
            water[st] = [("ideal device (L1)", lad["L1 ideal"][la]["F"]),
                         ("plain twin: IBM's gate, T1/T2 and readout record (L2)", lad["L2 plain twin (8a)"][la]["F"]),
                         ("realistic twin: + idle times, latency, asymmetric readout (L3 = 8c prediction)", P1["blocks"][k]["F"]["prediction"]),
                         ("L3 with the readout MEASURED in Round 1 (post-diction)", postdiction["variants"]["measured readout"][k]["mean"]),
                         ("hardware (Round 1)", r1["tomo"][k]["F"])]
        out["waterfall"] = water
    return out


def leg_T(r1: Mapping, kind: str) -> dict[int, np.ndarray]:
    """T matrices of the two legs from a Round 1 TOMO-9 result (keys int, or str after a JSON round trip)."""
    T = r1[kind]["T"]
    return {q: np.array(T[q] if q in T else T[str(q)]) for q in (2, 3)}


def _rz(a: float) -> np.ndarray:
    c, s = math.cos(a), math.sin(a)
    return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])


def dd_scenarios(r1: Mapping, r2: Mapping, cal: Mapping) -> dict:
    """What dynamical decoupling did during the 8 us wait (Part 9.8, P8.6 and Observation O8.7), as zero-parameter
    twirl-model predictions of the two DD blocks, each scored by chi^2 (2 values) against the hardware.

    T(N2) = Rc D: Rc = total rotation of the wait (+ native idle), D = incoherent remainder. IBM's DD refocuses Z
    phases but also pulses the idle message qubits, which are still in |0>; a pulsed neighbour in |0> adds the phase
    pi zeta t (zeta = static ZZ to the neighbour: Q41 on Q36, Q42 on Q43). Scenario S1 (DD keeps the non-Z part of the
    native rotation, spectators add their phase) fits; the no-spectator, same-sense and in-pair-ZZ readings do not."""
    from scipy.spatial.transform import Rotation
    from scipy.stats import chi2 as _chi2

    Tn, Tf = leg_T(r1, "calN2"), leg_T(r1, "fingerprint")
    Rc = {q: procrustes(Tn[q]).T for q in (2, 3)}
    Dn = {q: procrustes(Tn[q]) @ Tn[q] for q in (2, 3)}
    Rf = {q: procrustes(Tf[q]).T for q in (2, 3)}
    nonz = {}
    for q in (2, 3):
        v = Rotation.from_matrix(Rf[q]).as_rotvec(); v[2] = 0.0
        nonz[q] = Rotation.from_rotvec(v).as_matrix()
    ed = cal["edges"]
    s = {2: math.pi * ed["36-41"]["zz_kHz"] * 1e3 * TAU_US * 1e-6, 3: math.pi * ed["42-43"]["zz_kHz"] * 1e3 * TAU_US * 1e-6}
    zin = {2: 2 * math.pi * ed["21-36"]["zz_kHz"] * 1e3 * TAU_US * 1e-6, 3: 2 * math.pi * ed["43-44"]["zz_kHz"] * 1e3 * TAU_US * 1e-6}
    psi = message_state("product_HT")
    Fm = lambda M: leg_state_fidelity(psi, np.diag(np.diag(M[2])), np.diag(np.diag(M[3])))      # noqa: E731
    pair = lambda kept: {"DD": Fm({q: kept[q] @ Dn[q] for q in (2, 3)}),                       # noqa: E731
                         "TOMO+DD": Fm({q: Rc[q].T @ kept[q] @ Dn[q] for q in (2, 3)})}
    sc = {"S1: DD refocuses every Z phase, keeps the non-Z native rotation; pulsed spectators add pi zeta t":
          pair({q: _rz(s[q]) @ nonz[q] for q in (2, 3)}),
          "S0: as S1, but DD also removes the non-Z native rotation (the planning version)": pair({q: _rz(s[q]) for q in (2, 3)}),
          "S2: as S1, but DD leaves the whole native rotation": pair({q: _rz(s[q]) @ Rf[q] for q in (2, 3)}),
          "A1: no spectator phase (DD refocuses every Z phase only)": pair(nonz),
          "A0: no spectator phase, DD removes the whole rotation": pair({q: np.eye(3) for q in (2, 3)}),
          "S1': spectator phase in the same sense as the waiting rotation": pair({q: _rz(-s[q]) @ nonz[q] for q in (2, 3)}),
          "B: in-pair ZZ survives DD (the 8d reading)": pair({q: _rz(-zin[q]) @ nonz[q] for q in (2, 3)})}
    yd, sd = r2["tomo"]["R2_DDN2none_product_HT"]["F"], r2["tomo"]["R2_DDN2none_product_HT"]["F_sigma"]
    yt, st = r2["tomo"]["R2_DDN2corr_product_HT"]["F"], r2["tomo"]["R2_DDN2corr_product_HT"]["F_sigma"]
    rows = []
    for k, v in sc.items():
        c = ((v["DD"] - yd) / sd) ** 2 + ((v["TOMO+DD"] - yt) / st) ** 2
        rows.append({"scenario": k, "DD": v["DD"], "TOMO+DD": v["TOMO+DD"], "chi2": c, "p": float(_chi2.sf(c, 2))})
    return {"rows": rows, "spectator_rad": s, "hardware": {"DD": [yd, sd], "TOMO+DD": [yt, st]},
            "model_no_DD": {"none": Fm(Tn), "TOMO": Fm(Dn)}}


def n2_analysis(r1: Mapping, r2: Mapping, cal: Mapping) -> dict:
    """The 8 us wait seen by TOMO-9: T(N2) = Rz(phase) D with D the incoherent remainder. Scenarios for what dynamical
    decoupling (DD) did to the phase, a one-parameter fit (DD leaves a fraction g of each leg's phase) to the two DD
    blocks, and the incoherent part against the T1 / T2 record."""
    T, Tfp = leg_T(r1, "calN2"), leg_T(r1, "fingerprint")
    Dn = {q: procrustes(T[q]) @ T[q] for q in (2, 3)}
    Dfp = {q: procrustes(Tfp[q]) @ Tfp[q] for q in (2, 3)}
    th = {q: r1["calN2"][f"rotation_{q}"]["angle"] for q in (2, 3)}
    sgn = {q: float(np.sign(r1["calN2"][f"rotation_{q}"]["axis"][2])) for q in (2, 3)}

    def F_res(a2: float, a3: float, st: str = "product_HT") -> float:
        M = {q: np.diag(np.diag(_rz(a) @ Dn[q])) for q, a in ((2, a2), (3, a3))}
        return leg_state_fidelity(message_state(st), M[2], M[3])

    def Fs(key):
        t = (r1 if key.startswith("R1_") else r2)["tomo"][key]
        return t["F"], t.get("F_sigma") or 0.0

    hwF = {k: Fs(k) for k in ("R1_N2none_product_HT", "R2_N2corr_product_HT", "R2_DDN2none_product_HT", "R2_DDN2corr_product_HT",
                              "R1_N2none_bell_30_70", "R2_N2corr_bell_30_70")}
    ed = cal["edges"]
    zz = {2: 2 * math.pi * ed["21-36"]["zz_kHz"] * 1e3 * TAU_US * 1e-6, 3: 2 * math.pi * ed["43-44"]["zz_kHz"] * 1e3 * TAU_US * 1e-6}
    scen = {
        "no DD (measured channel)": {"DD": None, "TOMO+DD": None, "none": F_res(th[2], th[3]), "TOMO": F_res(0, 0)},
        "A: DD removes the whole waiting phase": {"DD": F_res(0, 0), "TOMO+DD": F_res(th[2], th[3])},
        "B: DD leaves the static ZZ phase, same sense as the total phase (8d interpretation)": {
            "DD": F_res(zz[2], zz[3]), "TOMO+DD": F_res(th[2] - zz[2], th[3] - zz[3])},
        "C: DD leaves the static ZZ phase, opposite sense": {"DD": F_res(zz[2], zz[3]), "TOMO+DD": F_res(th[2] + zz[2], th[3] + zz[3])},
    }
    yd, sd = hwF["R2_DDN2none_product_HT"]
    yt, st_ = hwF["R2_DDN2corr_product_HT"]
    gs = np.linspace(-0.6, 0.6, 1201)
    chi = [((F_res(g * th[2], g * th[3]) - yd) / sd) ** 2 + ((F_res((1 - g) * th[2], (1 - g) * th[3]) - yt) / st_) ** 2 for g in gs]
    g = float(gs[int(np.argmin(chi))])
    scen["D: fit, DD leaves a fraction g of the phase"] = {"g": g, "DD": F_res(g * th[2], g * th[3]),
                                                            "TOMO+DD": F_res((1 - g) * th[2], (1 - g) * th[3]),
                                                            "chi2": float(min(chi)), "leftover_rad": {2: g * th[2], 3: g * th[3]}}
    curve = [{"g": float(x), "DD": F_res(x * th[2], x * th[3]), "TOMO+DD": F_res((1 - x) * th[2], (1 - x) * th[3])} for x in gs[::50]]
    qrec = cal["qubits"]
    inc = {}
    for q, cq, bq in ((2, "q2", "q4"), (3, "q3", "q5")):
        inc[q] = {"measured_xyz": [Dn[q][a, a] / Dfp[q][a, a] for a in range(3)],
                  "record_xy": math.exp(-TAU_US / qrec[cq]["T2"]) * math.exp(-TAU_US / qrec[bq]["T2"]),
                  "record_z": math.exp(-TAU_US / qrec[cq]["T1"]) * math.exp(-TAU_US / qrec[bq]["T1"]),
                  "T1": (qrec[cq]["T1"], qrec[bq]["T1"]), "T2": (qrec[cq]["T2"], qrec[bq]["T2"])}
    k = lambda rad: rad / (2 * math.pi * TAU_US * 1e-6) / 1e3                 # noqa: E731  phase (rad) -> kHz
    return {"theta": th, "axis_z_sign": sgn, "kHz": {q: k(th[q]) for q in (2, 3)}, "zz_rad": zz,
            "zz_kHz": {2: ed["21-36"]["zz_kHz"], 3: ed["43-44"]["zz_kHz"]},
            "rest_kHz_same_sense": {q: k(th[q] - zz[q]) for q in (2, 3)},
            "rest_kHz_opposite_sense": {q: k(th[q] + zz[q]) for q in (2, 3)},
            "hardware": hwF, "scenarios": scen, "fit_curve": curve, "incoherent": inc,
            "D_N2": {q: Dn[q].tolist() for q in (2, 3)}, "D_fp": {q: Dfp[q].tolist() for q in (2, 3)},
            "model_tomo_bell": F_res(0, 0, "bell_30_70"), "model_none_bell": F_res(th[2], th[3], "bell_30_70")}
