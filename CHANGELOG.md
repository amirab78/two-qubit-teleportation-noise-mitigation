# Changelog

## 1.1.0 (2026-10-25)

The replication on real hardware: a new day and a second device. DOI:
[10.5281/zenodo.23639665](https://doi.org/10.5281/zenodo.23639665).

* Data: counts of eleven IBM jobs of 25 October 2026 on `ibm_kingston` and `ibm_fez`, with their calibration records,
  frozen circuits, pre-registered predictions (hashed before sending), approvals and the saved per-job analyses
  (`data/hardware/ibm_kingston_ibm_fez_2026-10-25/`, `data/reference/replication/`). Cleaned like the 4 October data
  (`scripts/clean_raw_data_2026_10_25.py`); the 4 October files are unchanged.
* `telecorr.hardware.replication`: the analysis of the eleven jobs from the raw counts: Round 1 / Round 2 on both
  devices, the waits at 4 / 8 / 16 µs, the DD spectator test, the same-job test on Fez, DD's rotation against the wait,
  and the correction learned under DD; model 2 (the twirl model) recomputed from the counts.
* `python -m telecorr.hardware reproduce --check` now covers both experiments: 8994 values, all equal to the saved
  results and to the pre-registered predictions. New option `--run`; `--out` writes one entry per experiment.
* Notebook 09 (replication); the roadmap of notebooks 01-08 links to it; notebook 08's open questions point to it.
* Tests: `tests/regression/test_replication.py` and the hash chain of the new data.
* Docs: `docs/replication.md` (the full record), the replication in the README, `docs/methods.md` and `data/README.md`.
* The citation now uses the concept DOI (all versions); v1.0.0 keeps its own DOI, 10.5281/zenodo.23634581.

## 1.0.0 (2026-10-25)

First public release. DOI: [10.5281/zenodo.23634581](https://doi.org/10.5281/zenodo.23634581).

* `telecorr` package: exact six-qubit protocol and qiskit circuits; channels and the Pauli twirl; probe cost and
  counts-only calibrators (COBYLA, automatic recalibration, Rotosolve); ceilings (Procrustes, two-qubit with
  certificate); process tomography (TOMO-9, TOMO-12, TOMO-75); crosstalk detector; Kalman tracker and the ABC rule;
  analysis of the IBM hardware runs from the raw counts.
* Notebooks 01-08: the study step by step, executed.
* Tests: the proven results of the study as unit tests; its saved numbers as regression tests (the hardware analysis
  to the last digit).
* Data: counts of the four IBM jobs on `ibm_kingston` (4 October 2026), the calibration record, the frozen circuits,
  the pre-registration record and the saved results.
