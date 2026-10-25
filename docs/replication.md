# The replication of 25 October 2026: full record

Eleven jobs on `ibm_kingston` (the device of 4 October) and `ibm_fez` (a second device, new chain), summarised in
notebook 9. Every value on this page is recomputed from the raw counts in `data/hardware/ibm_kingston_ibm_fez_2026-10-25/` by
`python -m telecorr.hardware reproduce --run 2026-10-25`; `--check` compares all of them (4283 values) with the
saved per-job analyses of the research project (`data/reference/replication/`) and with the pre-registered prediction
files, and `tests/regression/test_replication.py` does the same in the test suite.

## Procedure

The procedure of 4 October ([hardware.md](hardware.md)), with these differences:

1. **Two prediction models.** *Model 1* is the twin, a simulator built from IBM's calibration record with the
   measured waiting phase (not part of this repository: its values are shipped); it also gives each circuit's predicted
   distributions for the per-circuit check. *Model 2* is the Pauli-twirl model: each leg reaches Bob as diag(R T), with
   T the TOMO-9 matrix of the job's own Round 1. Model 2 is recomputed from the counts here and equals the values in
   the prediction files.
2. **Hashed before sending.** Each prediction file was written and its SHA-256 recorded (`preregistration/sha256.json`)
   before the jobs it covers were sent; each frozen manifest names its prediction file and that hash.
3. **Cost.** IBM's detailed estimate was taken on the live device for every job; the time cap was 1.25 × that estimate
   (at least 10 s). Every job was approved by its hash and cap before it was sent (`approvals.json`).
4. **Pilots first** on each device (K0, F0), and same-job designs for the follow-ups (G1, G2, G3), because the device
   drifts between jobs.

| job | purpose | device | IBM job id | circuits | shots | DD | IBM estimate (s) | approved cap (s) | QPU usage (s) | circuit time (s) | created |
|---|---|---|---|---|---|---|---|---|---|---|---|
| K0 | pilot | ibm_kingston | dbeqnffsstms73fciu80 | 12 | 6000 | no | 3.6 | 10 | 4 | 1.6 | 2026-10-25T07:11:57 UTC |
| K1 | Round 1 | ibm_kingston | dbeqouvgdj9s73fma3c0 | 130 | 349000 | no | 92.2 | 116 | 95 | 92.6 | 2026-10-25T07:15:07 UTC |
| K2 | Round 2 | ibm_kingston | dber02gm94ac738sb5tg | 54 | 126000 | no | 34.6 | 44 | 36 | 33.5 | 2026-10-25T07:30:19 UTC |
| K3 | DD, message qubits pulsed | ibm_kingston | dber2sfpc3tc739m7j00 | 27 | 63000 | yes | 18.6 | 24 | 19 | 17.0 | 2026-10-25T07:36:17 UTC |
| K4 | DD, message qubits not pulsed | ibm_kingston | dbere1nsstms73fck4pg | 27 | 63000 | yes | 18.6 | 24 | 19 | 17.0 | 2026-10-25T08:00:06 UTC |
| F0 | pilot | ibm_fez | dberne7gdj9s73fmbpd0 | 10 | 5000 | no | 3.3 | 10 | 4 | 1.3 | 2026-10-25T08:20:08 UTC |
| F1 | Round 1 | ibm_fez | dberor7pc3tc739m8tlg | 58 | 160000 | no | 42.6 | 54 | 44 | 41.4 | 2026-10-25T08:23:08 UTC |
| F2 | Round 2 | ibm_fez | dberrs8m94ac738scpig | 27 | 72000 | no | 20.3 | 27 | 21 | 18.6 | 2026-10-25T08:29:37 UTC |
| G1 | device alone and both corrections in one job | ibm_fez | dbes408m94ac738sd89g | 54 | 108000 | no | 29.4 | 38 | 30 | 28.0 | 2026-10-25T08:46:57 UTC |
| G2 | TOMO-9 under DD at 4 / 8 / 16 µs | ibm_kingston | dbes5lngdj9s73fmcl7g | 27 | 81000 | yes | 23.4 | 30 | 24 | 21.9 | 2026-10-25T08:50:30 UTC |
| G3 | DD + a correction learned under DD | ibm_kingston | dbes9k7pc3tc739m9vkg | 36 | 72000 | yes | 21.2 | 27 | 22 | 19.7 | 2026-10-25T08:58:56 UTC |

Total QPU time: 318 s.

## Pre-registered rules and criteria

| job | pre-registered rule or criterion | value | limit | outcome |
|---|---|---|---|---|
| K1 | go / no-go: per-circuit check: no circuit far from its prediction | 2 | 0 | **fails** |
| K1 | go / no-go: readout errors plausible | 0.043 | 0.080 | passes |
| K1 | go / no-go: native fidelity (product_HT) not far below the prediction | 0.891 | 0.814 | passes |
| K1 | go / no-go: native fidelity (bell_30_70) not far below the prediction | 0.857 | 0.783 | passes |
| K1 | go / no-go: TOMO-9 correction of N1 close to the exact inverse | 0.046 | 0.250 | passes |
| K1 | go / no-go: TOMO-9 predicts a large gain on N1 | 0.260 | 0.136 | passes |
| F1 | go / no-go: per-circuit check: no circuit far from its prediction | 0 | 0 | passes |
| F1 | go / no-go: readout errors plausible | 0.015 | 0.080 | passes |
| F1 | go / no-go: native fidelity (product_HT) not far below the prediction | 0.890 | 0.811 | passes |
| F1 | go / no-go: native fidelity (bell_30_70) not far below the prediction | 0.860 | 0.794 | passes |
| F1 | go / no-go: TOMO-9 correction of N1 close to the exact inverse | 0.142 | 0.250 | passes |
| F1 | go / no-go: TOMO-9 predicts a large gain on N1 | 0.276 | 0.132 | passes |
| K2 | REPLICATION (product_HT): F >= F(native) - 2 sigma | 0.889 | 0.873 | passes |
| K2 | REPLICATION (product_HT): recovered fraction >= 0.9 | 0.988 | 0.900 | passes |
| K2 | REPLICATION (bell_30_70): F >= F(native) - 2 sigma | 0.856 | 0.842 | passes |
| K2 | REPLICATION (bell_30_70): recovered fraction >= 0.9 | 0.998 | 0.900 | passes |
| K2 | WAIT 4us: F(corrected) - F(none) > 3 sigma | 6.044 | 3.000 | passes |
| K2 | WAIT 8us: F(corrected) - F(none) > 3 sigma | 19.713 | 3.000 | passes |
| K2 | WAIT 16us: F(corrected) - F(none) > 3 sigma | 48.931 | 3.000 | passes |
| F2 | REPLICATION (product_HT): F >= F(native) - 2 sigma | 0.910 | 0.872 | passes |
| F2 | REPLICATION (product_HT): recovered fraction >= 0.9 | 1.114 | 0.900 | passes |
| F2 | REPLICATION (bell_30_70): F >= F(native) - 2 sigma | 0.891 | 0.846 | passes |
| F2 | REPLICATION (bell_30_70): recovered fraction >= 0.9 | 1.089 | 0.900 | passes |
| G3 | G3-A 8us: F(DD + corr) - F(DD alone) > 3 sigma | -0.002 | σ = 0.011 | **fails** |
| G3 | G3-B 8us: DD + corr vs K2 fresh correction without DD (0.851) | -0.070 | σ = 0.011 | worse |
| G3 | G3-A 16us: F(DD + corr) - F(DD alone) > 3 sigma | 0.025 | σ = 0.012 | **fails** |
| G3 | G3-B 16us: DD + corr vs K2 fresh correction without DD (0.785) | -0.054 | σ = 0.011 | worse |

The one go / no-go failure (K1's per-circuit check) is explained under *Decisions disclosed*. The G3-A failures are a
result: a correction learned under DD barely helps.

## Every block

Round 1 blocks have no model-2 or twin column here: their predictions are in `x_pred_K0K1.json` and
`x_pred_F0F1_ibm_fez.json`. Bob's readout is mitigated with the matrix of the job's own Round 1 (K1 on Kingston, F1 on
Fez).

| job | block | message | hardware F ± σ | F, Bob's readout mitigated | model 2 | twin (model 1) |
|---|---|---|---|---|---|---|
| K1 | device alone | product_HT | 0.891 ± 0.006 | 0.898 | - | - |
| K1 | N1, no correction | product_HT | 0.720 ± 0.007 | 0.725 | - | - |
| K1 | N1, correction stored from 4 Oct | product_HT | 0.879 ± 0.006 | 0.885 | - | - |
| K1 | device alone | bell_30_70 | 0.857 ± 0.005 | 0.866 | - | - |
| K1 | N1, no correction | bell_30_70 | 0.532 ± 0.007 | 0.536 | - | - |
| K1 | N1, correction stored from 4 Oct | bell_30_70 | 0.865 ± 0.005 | 0.874 | - | - |
| K1 | 4 µs wait, no correction | product_HT | 0.808 ± 0.008 | 0.814 | - | - |
| K1 | 8 µs wait, no correction | product_HT | 0.626 ± 0.009 | 0.630 | - | - |
| K1 | 16 µs wait, no correction | product_HT | 0.261 ± 0.007 | 0.261 | - | - |
| F1 | device alone | product_HT | 0.890 ± 0.006 | 0.898 | - | - |
| F1 | N1, no correction | product_HT | 0.706 ± 0.007 | 0.711 | - | - |
| F1 | device alone | bell_30_70 | 0.860 ± 0.005 | 0.871 | - | - |
| F1 | N1, no correction | bell_30_70 | 0.511 ± 0.007 | 0.516 | - | - |
| K2 | N1, TOMO-9 correction | product_HT | 0.889 ± 0.006 | 0.896 | 0.879 | 0.948 |
| K2 | N1, TOMO-9 correction | bell_30_70 | 0.856 ± 0.005 | 0.866 | 0.850 | 0.917 |
| K2 | 4 µs wait, TOMO-9 correction | product_HT | 0.877 ± 0.008 | 0.884 | 0.862 | 0.914 |
| K2 | 8 µs wait, TOMO-9 correction | product_HT | 0.851 ± 0.007 | 0.857 | 0.825 | 0.874 |
| K2 | 16 µs wait, TOMO-9 correction | product_HT | 0.785 ± 0.008 | 0.791 | 0.797 | 0.825 |
| F2 | N1, TOMO-9 correction | product_HT | 0.910 ± 0.006 | 0.919 | 0.899 | 0.953 |
| F2 | N1, TOMO-9 correction | bell_30_70 | 0.891 ± 0.005 | 0.903 | 0.863 | 0.941 |
| K3 | DD alone, message pulsed (8 µs) | product_HT | 0.778 ± 0.008 | 0.783 | 0.820 | 0.880 |
| K3 | K1 correction + DD, message pulsed (8 µs) | product_HT | 0.524 ± 0.008 | 0.526 | 0.580 | 0.667 |
| K4 | DD alone (8 µs) | product_HT | 0.815 ± 0.008 | 0.821 | 0.825 | 0.878 |
| K4 | K1 correction + DD (8 µs) | product_HT | 0.583 ± 0.008 | 0.586 | 0.635 | 0.672 |
| G1 | device alone | product_HT | 0.907 ± 0.008 | 0.915 | 0.901 | 0.954 |
| G1 | N1, TOMO-9 correction | product_HT | 0.909 ± 0.008 | 0.918 | 0.899 | 0.950 |
| G1 | device alone, fingerprint correction | product_HT | 0.908 ± 0.008 | 0.916 | 0.906 | 0.950 |
| G1 | device alone | bell_30_70 | 0.857 ± 0.006 | 0.868 | 0.875 | 0.943 |
| G1 | N1, TOMO-9 correction | bell_30_70 | 0.880 ± 0.006 | 0.892 | 0.863 | 0.929 |
| G1 | device alone, fingerprint correction | bell_30_70 | 0.896 ± 0.006 | 0.908 | 0.880 | 0.939 |
| G3 | DD alone (8 µs) | product_HT | 0.783 ± 0.008 | 0.788 | 0.785 | 0.883 |
| G3 | correction learned under DD + DD (8 µs) | product_HT | 0.780 ± 0.008 | 0.786 | 0.788 | 0.874 |
| G3 | DD alone (16 µs) | product_HT | 0.705 ± 0.009 | 0.710 | 0.706 | 0.813 |
| G3 | correction learned under DD + DD (16 µs) | product_HT | 0.731 ± 0.008 | 0.736 | 0.719 | 0.788 |

## The waiting phase (K1)

The TOMO-9 rotation at 4, 8 and 16 µs is about −z on both legs. A fit through the origin, θ = 2π f t, gives
f = 15.78 kHz (leg q2) and 15.10 kHz (leg q3), with residuals below
0.026 rad. On 4 October the 8 µs wait gave
12.77 / 22.33 kHz: the frequency drifts from day to day.

## DD: the spectator test (K3, K4) and the rotation against the wait (G2)

Model 2 for K3 / K4 is Part 9.8's: DD refocuses every z phase and keeps the non-z native rotation; pulsed spectators
add R_z(π ζ t), with ζ from the K1 calibration record.

| job | leg | z rotation under DD (rad) | model 2 | pull |
|---|---|---|---|---|
| K3 | q2 | 0.199 ± 0.016 | 0.114 | +5.2 |
| K3 | q3 | 0.191 ± 0.017 | 0.127 | +3.7 |
| K4 | q2 | 0.091 ± 0.016 | 0.000 | +5.8 |
| K4 | q3 | 0.081 ± 0.016 | 0.000 | +5.0 |

G2: z(t) = a + b t on each leg (CONST if b is consistent with 0, LINEAR if a is, each at 2σ):

| leg | intercept a (rad) | slope b (rad / µs) | b / σ | χ² (1 dof) | reading |
|---|---|---|---|---|---|
| q2 | -0.027 ± 0.019 | 0.016 ± 0.002 | 8.237 | 0.989 | LINEAR |
| q3 | -0.065 ± 0.020 | 0.017 ± 0.002 | 8.165 | 0.718 | undecided |

## Same-job test on Fez (G1)

REAL = the two-job gap F2 − F1; DRIFT = model 2's gap (both written down before G1 was sent).

| message | same-job difference | D ± σ | D / σ | REAL | DRIFT | log10 likelihood ratio REAL / DRIFT |
|---|---|---|---|---|---|---|
| product_HT | N1 + correction - native (same job) | 0.002 ± 0.011 | 0.2 | 0.021 | -0.002 | -0.559 |
| product_HT | fingerprint correction - native (same job) | 0.001 ± 0.011 | 0.1 | - | - | - |
| bell_30_70 | N1 + correction - native (same job) | 0.024 ± 0.008 | 2.8 | 0.031 | -0.012 | 3.710 |
| bell_30_70 | fingerprint correction - native (same job) | 0.039 ± 0.009 | 4.5 | - | - | - |

## Per-circuit checks

Every circuit against the twin's predicted distributions, with the rule of 4 October. The twin has no DD offset, so it
was expected to flag DD circuits.

| job | circuits | flagged | flagged circuits |
|---|---|---|---|
| K0 | 12 | 0 | - |
| K1 | 130 | 2 | `R1_N2none16us_product_HT_ZX`, `R1_N2none16us_product_HT_ZY` |
| K2 | 54 | 2 | `R2_ver_tomo9_s0_X`, `R2_ver_tomo9_s0_Y` |
| K3 | 27 | 9 | `R3_DDcal_tomo9_s0_Z`, `R3_DDcal_tomo9_s0_X`, `R3_DDcal_tomo9_s0_Y`, `R3_DDcal_tomo9_s2_X`, `R3_DDnone_product_HT_XX`, `R3_DDcorr_product_HT_XX`, `R3_DDcorr_product_HT_XY`, `R3_DDcorr_product_HT_YX`, `R3_DDcorr_product_HT_YY` |
| K4 | 27 | 2 | `R4_DDUcal_tomo9_s0_X`, `R4_DDUcal_tomo9_s0_Y` |
| F0 | 10 | 0 | - |
| F1 | 58 | 0 | - |
| F2 | 27 | 0 | - |
| G1 | 54 | 0 | - |
| G2 | 27 | 7 | `R6_DDU4uscal_tomo9_s0_X`, `R6_DDU8uscal_tomo9_s0_X`, `R6_DDU8uscal_tomo9_s0_Y`, `R6_DDU16uscal_tomo9_s0_Z`, `R6_DDU16uscal_tomo9_s0_X`, `R6_DDU16uscal_tomo9_s0_Y`, `R6_DDU16uscal_tomo9_s2_X` |
| G3 | 36 | 2 | `R7_DDU8usnone_product_HT_YX`, `R7_DDU16usnone_product_HT_YX` |

## Crosstalk detector

| job | TOMO-9 block | X² | dof | p |
|---|---|---|---|---|
| K1 | cal9 | 23.933 | 9 | 0.004 |
| K1 | fp | 55.845 | 8 | 3.02e-09 |
| K1 | calN24us | 15.721 | 9 | 0.073 |
| K1 | calN28us | 11.961 | 9 | 0.215 |
| K1 | calN216us | 27.149 | 9 | 0.001 |
| F1 | cal9 | 10.621 | 9 | 0.303 |
| F1 | fp | 4.140 | 7 | 0.763 |
| K2 | ver | 39.940 | 8 | 3.29e-06 |
| F2 | ver | 8.244 | 6 | 0.221 |

The leg-local TOMO-9 model is rejected at the fingerprint and at the largest signals, as on 4 October: the chain has a
small correlated error between the legs. It cost nothing measurable in the corrected fidelities.

## Decisions disclosed

* **K1's per-circuit check flagged two circuits** (both at 16 µs). As pre-registered, the run stopped and the circuits
  were inspected: the waiting frequency had drifted since 4 October, and the twin with K1's fitted frequency flags none
  (`data/reference/replication/x_K1_retro_wait.json`). The run continued, and **K2's predictions use the frequency K1
  measured**: a decision made after seeing data, before K2 was sent.
* **F2's time cap** got a floor of 1.2 × the measured-rate estimate (26 → 27 s; the circuits are unchanged).
* **G1's ideal-simulator test band** was widened from ±0.02 to ±0.03, which matches the shot noise at 2000 shots.
* None of these changed a prediction or a success criterion after its job was sent.
* **Drift between jobs.** DD alone at 8 µs gave 0.815 in K4 and 0.783 in G3, the same circuits:
  comparisons between jobs carry about ±0.03 beyond the shot-noise σ, which is why G1, G2 and G3 are same-job designs.

## Not reproduced here

* The twin's predictions (model 1), the per-circuit predicted distributions and the K1 retrodiction are shipped as
  saved values; the twin is not part of this repository.
* The tests on an ideal simulator run before each send (`x_l1_tests_*.json`) and the cost estimates on the live device
  are shipped as saved values.
