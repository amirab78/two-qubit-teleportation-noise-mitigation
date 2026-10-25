# Data

Everything needed to rebuild the hardware results of this repository, and the saved results the notebooks and tests
compare with. 48 files, 2.3 MB in total; every file is listed in `MANIFEST.json` with its source in the
research project, the SHA-256 of the source and of the shipped file, and the transformation applied.

| folder | files | size |
|---|---|---|
| `data/hardware/ibm_kingston_2026-10-04/calibration/` | 1 | 0.06 MB |
| `data/hardware/ibm_kingston_2026-10-04/circuits/` | 12 | 0.46 MB |
| `data/hardware/ibm_kingston_2026-10-04/counts/` | 4 | 0.73 MB |
| `data/hardware/ibm_kingston_2026-10-04/preregistration/` | 11 | 0.45 MB |
| `data/reference/` | 14 | 0.56 MB |
| `data/reference/simulation/` | 6 | 0.05 MB |

## The hardware experiment (`hardware/ibm_kingston_2026-10-04/`)

Four jobs on the IBM device `ibm_kingston` on 4 October 2026 (SamplerV2, no dynamical decoupling or twirling unless
stated, every option set explicitly):

| job | IBM job id | shots | QPU usage (s) | approved cap (s) | circuit time (s) | queue (s) | created |
|---|---|---|---|---|---|---|---|
| pilot (R0) | db16kn2vog1s73fhnpog | 5000 | 4 | 10 | 1.3 | 54 | 2026-10-04T15:01:48 UTC |
| Round 1 | db16p1bid5ic73eqnsug | 436000 | 117 | 201 | 114.9 | 99 | 2026-10-04T15:11:01 UTC |
| Round 2 main | db175bjid5ic73eqoav0 | 210000 | 58 | 99 | 55.5 | 143 | 2026-10-04T15:37:18 UTC |
| Round 2 DD | db175givog1s73fhocjg | 54000 | 17 | 28 | 14.6 | 125 | 2026-10-04T15:37:38 UTC |

* `counts/<job id>.json.gz`: one file per job. `pubs` is a list of PUB results in the order of the frozen job, each
  with its `label`, the `counts` per classical register (`alice`: 4 bits, Alice's measurement; `bob`: 2 bits, keys
  `"b5b4"`; `msg`: the message measured directly) and `shots_int`, the outcome of every shot as an integer (needed
  for the Pauli frame and the bootstrap). PUBs with a parameter array hold one entry per binding. `record` holds
  the job's experimental fields (device, shots, time cap, hash, submission time, usage); `metrics` holds IBM's time
  stamps and circuit execution time.
* `calibration/`: the device calibration record (T1, T2, readout, gate errors, static ZZ) in force when the jobs ran;
  the four job-time records were identical, one is shipped.
* `circuits/`: the frozen jobs, transpiled for the device: QPY circuits, parameter values of every binding, and a
  manifest (labels, binding index, shots, time cap, hash).
* `preregistration/`: the twin's predictions and their SHA-256 records, the Round 2 addendum and twirl-model prediction,
  the configuration choices, the approvals (hash and time cap of each job) and the job ledger.

**Circuit labels.** A binding label is `<block>_<message>_<bases>`, e.g. `R2_N1corr_bell_30_70_XY` (Bob measured in X on
q4 and Y on q5). Other labels: `readout_<pattern>` (readout calibration), `alice_<message>_<bases>` (the message measured
directly), `..._unprep` (un-preparation check), `tomo9_s<j>_<axis>` and `tomo12_...` (TOMO probe circuits, setting j,
probe axis), `R2_ver_...` (the Round 2 TOMO-9 verification). The block keys:

| block key | meaning (version A unless stated) |
|---|---|
| `R1_nativeA` | device alone |
| `R1_nativeB` | device alone, version B |
| `R1_C` | device alone, version C |
| `R1_N1none` | N1, no correction |
| `R1_N1noneB` | N1, no correction, version B |
| `R1_N1oracle` | N1, exact inverse (oracle) |
| `R2_N1corr` | N1, TOMO-9 correction |
| `R2_N1corrB` | N1, TOMO-9 correction, version B |
| `R2_N0corr` | device alone, TOMO-9 correction |
| `R1_N2none` | 8 µs wait, no correction |
| `R2_N2corr` | 8 µs wait, TOMO-9 correction |
| `R2_DDN2none` | 8 µs wait, DD |
| `R2_DDN2corr` | 8 µs wait, TOMO-9 correction + DD |

Versions: A = Pauli frame (all measurements at the end, Bob's corrections in software); B = mid-circuit measurement
and feed-forward; C = deferred measurement.

## Saved results (`reference/`)

The research project's saved results that the regression tests and notebooks compare with: the hardware analysis
(`8d_*`, `8e_*`, `9_8_*`), the twin's ladder and post-diction (`8c_ladder.json`, `8d_postdiction.json`), the project's
table of final numbers (`final_numbers.json`) and, in `simulation/`, the saved simulation results used by notebooks 3,
6 and 7.

## What was changed

The counts and per-shot outcomes are unchanged (checked pub by pub when the data were made). Removed: every field
naming the account, instance, plan quota or billing, and the id of a usage-detail metric. Converted: the time stamps
of the calibration record, from the local time zone of the recording computer to UTC. Everything else is a byte copy.
`scripts/clean_raw_data.py` documents how this folder was made from the research project.

## Use and citation

The data are released with the code under the MIT licence. Please cite the repository (see `CITATION.cff`).
