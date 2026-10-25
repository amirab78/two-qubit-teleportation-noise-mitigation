# The hardware experiment: full record

This page holds the details of the IBM experiment of 4 October 2026 that notebook 8 summarises: the jobs, every
pre-registered hypothesis, every block against its prediction, and the per-circuit checks. Every value is recomputed
from the raw counts in `data/` by `python -m telecorr.hardware reproduce --run 2026-10-04`; `--check` compares all of
them (4711 values) with the saved results of the research project, and `tests/regression/test_hardware.py` does
the same in the test suite. The replication of 25 October (notebook 9) has its own record:
[replication.md](replication.md).

## Procedure

1. **Freeze.** Each job was built once, saved (QPY circuits, parameter arrays, a manifest with labels, shots and a
   SHA-256 hash of the whole job) and never rebuilt (`data/hardware/ibm_kingston_2026-10-04/circuits/`).
2. **Predict.** The outcome of every block was predicted by a simulator built from IBM's calibration record (the
   *twin*; not part of this repository), with a band, and the prediction files were hashed before sending
   (`data/hardware/ibm_kingston_2026-10-04/preregistration/`). The Round 2 corrections were added in a hashed addendum, and the twirl-model
   prediction of Round 2 was hashed before Round 2 was sent.
3. **Approve.** Each job was approved by its hash, with an exact time cap (`approvals.json`). The send guard refused
   any job that was not approved, had another cap, or had been sent before (`scripts/submit_hardware_jobs.py`
   documents the procedure).
4. **Run and analyse.** Four jobs on `ibm_kingston`, 4 October 2026:

| job | IBM job id | shots | QPU usage (s) | approved cap (s) | circuit time (s) | queue (s) | created |
|---|---|---|---|---|---|---|---|
| pilot (R0) | db16kn2vog1s73fhnpog | 5000 | 4 | 10 | 1.3 | 54 | 2026-10-04T15:01:48 UTC |
| Round 1 | db16p1bid5ic73eqnsug | 436000 | 117 | 201 | 114.9 | 99 | 2026-10-04T15:11:01 UTC |
| Round 2 main | db175bjid5ic73eqoav0 | 210000 | 58 | 99 | 55.5 | 143 | 2026-10-04T15:37:18 UTC |
| Round 2 DD | db175givog1s73fhocjg | 54000 | 17 | 28 | 14.6 | 125 | 2026-10-04T15:37:38 UTC |

`tests/regression/test_data_integrity.py` checks the hash chain: approved hash = frozen manifest hash = the hash in
the job record, and every prediction file was hashed before its job was created.

## Pre-registered hypotheses

|  | hypothesis as pre-registered | evidence (hardware) | verdict |
|---|---|---|---|
| H1 | native F inside its band; fingerprint below the 0.05 rad floor (no N0 block) | native A 0.875 / 0.861 (inside: False / False); fingerprint 0.159 rad | fails |
| H2 | TOMO-9 from real counts recovers the whole N1 loss (twin 1.00 +/- 0.06 / 1.00 +/- 0.02) | product_HT: 1.057 +/- 0.057; bell_30_70: 0.975 +/- 0.020 | holds |
| H3 | TOMO-9 = oracle within noise (TOMO-9 > oracle only with native coherent errors on the channel) | product_HT: TOMO-9 - oracle = +0.0041 +/- 0.0092; bell_30_70: TOMO-9 - oracle = +0.0029 +/- 0.0069 | holds |
| H4 | N2: only the coherent part is recoverable; TOMO-9 alone ~ DD alone; TOMO-9 + DD over-corrects | none 0.562, tomo9 0.809, dd 0.775, tomo9_dd 0.499 | descriptive |
| H5/H10 | readout biases the fidelity, not the correction (A); mitigated - raw ~ +0.01 in A | product_HT: mitigated - raw = +0.0088; bell_30_70: mitigated - raw = +0.0128 | holds |
| H6 | Bob-mitigated A > C > B, all within ~0.01 (product_HT) | A 0.884, C 0.919, B 0.882 | fails |
| H7 | the twin is optimistic | 12 of 12 Round 1 blocks below the L3 prediction | holds |
| H8 | the crosstalk detector is silent | p = 0.00021 (fingerprint), 0.845 (N1 calibration) | fails |
| H9 | a native coherent rotation is claimed only above ~0.05 rad | fingerprint 0.159 rad | descriptive |
| H11 | the A - B gap measures the feed-forward wait (twin: 0.011 / 0.011 mitigated) | product_HT: A - B (mitigated) = +0.0020; bell_30_70: A - B (mitigated) = +0.0111 | descriptive |

"Descriptive" hypotheses had no pass / fail rule. The twin was optimistic in absolute level (H7); the relative claims
about the correction (H2, H3, H5/H10) held; the crosstalk detector fired on the device alone (H8), which is an open
question.

## Every block against its prediction

| round | block (version A unless stated) | message | hardware F ± σ | predicted F (twin) | predicted band | inside |
|---|---|---|---|---|---|---|
| 1 | device alone | product_HT | 0.8751 ± 0.0063 | 0.9674 | [0.911, 0.990] | no |
| 1 | device alone | bell_30_70 | 0.8608 ± 0.0046 | 0.9503 | [0.907, 0.963] | no |
| 1 | device alone, version B | product_HT | 0.8819 ± 0.0068 | 0.9660 | [0.927, 0.984] | no |
| 1 | device alone, version B | bell_30_70 | 0.8460 ± 0.0047 | 0.9370 | [0.891, 0.952] | no |
| 1 | N1, no correction | product_HT | 0.7109 ± 0.0065 | 0.7646 | [0.718, 0.790] | no |
| 1 | N1, no correction | bell_30_70 | 0.5217 ± 0.0069 | 0.5732 | [0.528, 0.609] | no |
| 1 | N1, exact inverse (oracle) | product_HT | 0.8804 ± 0.0060 | 0.9605 | [0.911, 0.988] | no |
| 1 | N1, exact inverse (oracle) | bell_30_70 | 0.8495 ± 0.0049 | 0.9503 | [0.906, 0.961] | no |
| 1 | 8 µs wait, no correction | product_HT | 0.5623 ± 0.0062 | 0.8977 | [0.680, 0.918] | no |
| 1 | 8 µs wait, no correction | bell_30_70 | 0.5935 ± 0.0069 | 0.8641 | [0.684, 0.878] | no |
| 1 | N1, no correction, version B | bell_30_70 | 0.5193 ± 0.0071 | 0.5781 | [0.524, 0.602] | no |
| 1 | device alone, version C | product_HT | 0.9185 ± 0.0072 | 0.9720 | [0.930, 0.999] | no |
| 2 | N1, TOMO-9 correction | product_HT | 0.8845 ± 0.0070 | 0.9633 | [0.913, 0.983] | no |
| 2 | N1, TOMO-9 correction | bell_30_70 | 0.8523 ± 0.0049 | 0.9479 | [0.906, 0.967] | no |
| 2 | 8 µs wait, TOMO-9 correction | product_HT | 0.8089 ± 0.0062 | 0.8966 | [0.754, 0.923] | yes |
| 2 | 8 µs wait, TOMO-9 correction | bell_30_70 | 0.7643 ± 0.0059 | 0.8613 | [0.742, 0.879] | yes |
| 2 | N1, TOMO-9 correction, version B | bell_30_70 | 0.8519 ± 0.0050 | 0.9370 | [0.894, 0.952] | no |
| 2 | device alone, TOMO-9 correction | product_HT | 0.8726 ± 0.0062 | - | - | - |
| 2 | device alone, TOMO-9 correction | bell_30_70 | 0.8556 ± 0.0049 | - | - | - |
| 2 | 8 µs wait, DD | product_HT | 0.7749 ± 0.0070 | 0.8955 | [0.749, 0.916] | yes |
| 2 | 8 µs wait, TOMO-9 correction + DD | product_HT | 0.4994 ± 0.0067 | 0.8964 | [0.680, 0.916] | no |

## Per-circuit checks

Every circuit's measured distributions (Alice's 16 outcomes and Bob's 4, Bob frame-corrected) were compared with the
pre-registered prediction by the total variation distance. A circuit is flagged when the distance exceeds the 99.9 %
quantile of the shot noise at its number of shots plus max(0.05, the prediction's band).

| round | circuits | within their threshold | flagged |
|---|---|---|---|
| 1 | 172 | 149 | readout_101010, R1_nativeB_unprep, R1_N1oracle_unprep, R1_nativeB_bell_30_70_YY, R1_N1noneB_bell_30_70_YY, R1_fp_tomo9_s0_X, R1_fp_tomo9_s0_Y, R1_cal9_tomo9_s0_Y, R1_cal12_tomo12_s2_Y, R1_N2none_product_HT_XX, R1_N2none_product_HT_XY, R1_N2none_product_HT_XZ, R1_N2none_product_HT_YX, R1_N2none_product_HT_YY, R1_N2none_product_HT_YZ, R1_N2none_product_HT_ZX, R1_N2none_product_HT_ZY, R1_N2none_bell_30_70_XX, R1_N2none_bell_30_70_YY, R1_calN2_tomo9_s0_X, R1_calN2_tomo9_s0_Y, R1_calN2_tomo9_s1_Y, R1_calN2_tomo9_s2_X |
| 2 | 91 | 83 | R2_N2corr_product_HT_XX, R2_N2corr_product_HT_XY, R2_N2corr_product_HT_YX, R2_N2corr_product_HT_YY, R2_N2corr_product_HT_ZX, R2_N2corr_product_HT_ZY, R2_N2corr_bell_30_70_XX, R2_N2corr_bell_30_70_YY |

## Not reproduced here

* The go / no-go rules applied between the rounds used the remaining account quota at that moment (account data,
  not shipped).
* The twin's predictions and the post-diction with the measured readout are shipped as saved results
  (`data/reference/`); the twin itself (a device simulator built from IBM's calibration record) is not part of this
  repository.
* IBM's own time estimate of each job needs the device's transpile target at that time; the estimates recorded while
  the jobs were queued are in the job records.
