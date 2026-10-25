# Changelog

## 1.0.0 (2026-10-25)

First public release.

* `telecorr` package: exact six-qubit protocol and qiskit circuits; channels and the Pauli twirl; probe cost and
  counts-only calibrators (COBYLA, automatic recalibration, Rotosolve); ceilings (Procrustes, two-qubit with
  certificate); process tomography (TOMO-9, TOMO-12, TOMO-75); crosstalk detector; Kalman tracker and the ABC rule;
  analysis of the IBM hardware runs from the raw counts.
* Notebooks 01-08: the study step by step, executed.
* Tests: the proven results of the study as unit tests; its saved numbers as regression tests (the hardware analysis
  to the last digit).
* Data: counts of the four IBM jobs on `ibm_kingston` (4 October 2026), the calibration record, the frozen circuits,
  the pre-registration record and the saved results.
