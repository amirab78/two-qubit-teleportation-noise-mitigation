# Methods, notation and limitations

A compact reference for the methods implemented in `telecorr`, in the order of the eight notebooks. The proofs are in
the theory reports of the research project. The claim labels (T = theorem, P = proposition, O = observation) follow
them: T5.6 is result 6 of Part 9.5, which belongs to step 5. `docs/claims.md` lists the test that checks each one.

## Notation

| symbol | meaning |
|---|---|
| q0, q1 | Alice's message qubits |
| q2, q3 | channel qubits (Alice's halves of the Bell pairs (q2, q4) and (q3, q5)) |
| q4, q5 | Bob's qubits; Bob's outcome keys are written `b5b4` (qiskit order) |
| `product_HT`, `bell_30_70` | the two test messages: H then T on each qubit; √0.3\|00⟩ + √0.7\|11⟩ |
| N1 | the coherent channel noise: rotations by 0.5 rad (q2) and 1.3 rad (q3) |
| T | the 3 × 3 transfer matrix of one channel leg (unital block of its Pauli transfer matrix); 15 × 15 for two qubits |
| R, R_V | the Bloch rotation of a correction V |
| λ = diag(R T) | the shrink factors of the Pauli channel that reaches Bob through one leg |
| p_I = (1 + tr(R T))/4 | probability that a leg delivers its qubit unchanged |
| F, F_e, F_avg | fidelity to one message; entanglement fidelity; Haar average over all two-qubit messages, F_avg = (4 F_e + 1)/5 |
| ceiling | the largest F_avg any correction of the given form can reach |
| gap | ceiling − F_avg |
| P(00) | the un-preparation check: Bob undoes the message preparation and reads 00; equals F for a pure message |
| C | the probe cost, 1 − mean of P(00) over the probes \|00⟩, \|++⟩, \|+i+i⟩ |
| TOMO-9 / TOMO-12 / TOMO-75 | process tomography from counts with 9 / 12 circuits per leg, or 75 two-qubit circuits |
| N0, N1, N2 | hardware blocks: the device alone, injected rotations, an 8 µs wait |
| versions A, B, C | Pauli frame in software; mid-circuit measurement with feed-forward; deferred measurement |
| twin | a simulator of the device built from IBM's calibration record (not part of this repository) |
| DD | dynamical decoupling |

## The methods, step by step

**1. Teleportation** (`protocol`). The protocol is simulated exactly as a 64 × 64 density matrix, branch by branch
over Alice's 16 outcomes, with Bob's Pauli corrections in each branch. With a perfect channel it is the identity
channel (T1.5), and Alice's outcomes are uniform and independent of the message (T1.6). The same protocol is available
as qiskit circuits (deferred or measured) for Aer.

**2. Channel noise** (`channels`, `twirl`). Channels are given by Kraus operators, Choi matrices or transfer matrices.
Averaging over Alice's outcomes twirls whatever acts on a channel qubit into the Pauli channel with shrink factors
diag(T) (T2.4, T5.2). Every branch is pure, but Bob's average state is mixed (T2.5). After teleportation no
operation on Bob's side helps an unknown message (T2.14). A correction must act on the channel qubits.

**3. Correction from counts** (`probes`, `calibrate.cobyla`). Six angles U(θ, φ, λ) on q2, q3 are chosen by COBYLA to
minimise the measured probe cost C. The probe probabilities are products of the leg factors (1 + λ_a)/2 (T3.5), and
(5/6)(1 − F_avg) ≤ C ≤ (5/6)(1 − F_avg) + (25/48)(1 − F_avg)² (T3.6, T5.9). One probe is blind to rotations about
its own axis (T3.8). Counts are drawn from the exact multinomial law of the circuits (`ProbeSampler`), which equals
the Aer distribution of the qiskit circuits.

**4. Automatic recalibration** (`calibrate.AutoCalibrator`). Monitor the probe P(00), recalibrate below a threshold
with a multi-start search, verify. The cost is a bilinear form in the squared quaternions of the two residual
rotations (P4.1). It has no spurious local minimum (T4.2), and its other critical points are saddles at
C ∈ {2/3, 5/6, 8/9, 1} (T4.3). Single starts fail through the poorly conditioned U3 chart, the evaluation budget,
shot noise and rare saddle traps. An accepted start with C < 0.01 guarantees F > 0.985 for every message (P4.10).

**5. Incoherent noise and the ceiling** (`ceiling`, `calibrate.rotosolve`). The best correction maximises tr(R T)
over rotations: orthogonal Procrustes, closed form from the SVD with a determinant fix. It is the inverse of the
rotation exactly when the incoherent part is symmetric and the sum of any two of its eigenvalues is non-negative (T5.6). Rotosolve updates one local
coordinate at a time with an exact sinusoid fit (P5.15). Under non-depolarizing noise the probe cost's minimum is
slightly off the ceiling. This bias is real but second order (P5.10).

**6. Correlated noise** (`correlated`, `detector`). The two-qubit twirl gives λ = diag(Ad(W) T) (T6.2) and
F_avg = (4 p_II + 1)/5 (T6.3). The best local correction comes from alternating Procrustes. The best two-qubit
correction comes from a χ-matrix start and Newton steps on SU(4), certified by p_II ≤ λ_max(χ). The inverse of the
unitary part is optimal when the incoherent part is symmetric positive semidefinite (T6.5); otherwise it is certified
channel by channel. TOMO-75 reads every entry of the 15 × 15 matrix exactly once (T6.11). Crosstalk is detected for
free from TOMO-9 data by a Pearson test of independence of Bob's two bits (P6.13).

**7. Efficient learning and tracking** (`tomography`, `tracking`). TOMO-9 reads every entry of each leg's T once
(T7.2), and Procrustes of the estimate gives the correction (T7.3). The gap falls as c / shots with c ≈ 4.9 (P7.5). A
Kalman filter on vec(T) tracks a drifting channel (P7.8). The ABC rule spends the fewest shots that keep the predicted
gap below a target. Its saving over repeated TOMO-9 follows from information balance (P7.9).

**8. Real hardware** (`hardware`). Version A applies Bob's corrections in software: flip(Z) = x, flip(X) = z,
flip(Y) = x ⊕ z, which equals correcting before the measurement (T8.1). The readout law is λ̂ = (1 − e0 − e1) λ + (e1 − e0)
(T8.2). Bob's state is estimated by linear tomography (unbiased, exact multinomial variance; T8.3) with a parametric
bootstrap. The recovered fraction f = (F_corr − F_none)/(F_native − F_none) has a delta-method error and an exact Fieller
interval (P8.10).

## Statements corrected during the study

The repository uses the corrected versions only.

* Single-start failures are **not** caused by local minima: the probe cost has none (T4.2). They come from the U3 angle
  chart, the evaluation budget, shot noise and a rare saddle trap.
* The probe bias under non-depolarizing noise is real but second order (P5.10), and zero for depolarizing noise.
* "The inverse of the unitary part is the best two-qubit correction" is a theorem only for symmetric positive
  semidefinite incoherent parts (T6.5). A counterexample exists for non-unital noise. For the study's channels it is
  certified by a bound.
* The large-N formula for the crosstalk detector's power is wrong at small N, because the skip rule drops the most
  informative table.
* Two rounds of tomography cost 15.44 / shots: re-centring between rounds **hurts**.
* The dynamical-decoupling shortfall on hardware is not caused by the in-pair ZZ phase. IBM's DD also pulses the idle
  message qubits, and a pulsed neighbour in |0⟩ adds a phase through the static ZZ coupling (P8.6, O8.7). This
  explanation has no free parameter.
* A doubly damped Bell pair keeps ⟨ZZ⟩ = (1 − γa)(1 − γb) + γa γb.
* For correlated noise, the values 0.8341 / 0.8291 (best local) and 0.9112 / 0.9073 (two-qubit) are **ceilings**. The
  corrections learned from counts in the study reached 0.8338 / 0.8288 (local TOMO-9) and 0.9112 / 0.9072 (TOMO-75).

## Limitations

* **One device, one day.** The hardware experiment ran on `ibm_kingston` on 4 October 2026 (196 s of QPU time). It
  should be repeated on other days and devices.
* **The twin was optimistic.** The simulator built from IBM's calibration record missed the absolute hardware
  fidelities. The relative predictions and the TOMO-9 twirl model held.
* **Not run on hardware:** the two-qubit correction (TOMO-75) and drift tracking (Kalman / ABC).
* **Simulation assumptions.** Steps 2–7 use model noise families, with ideal state preparation and measurement in
  simulation. Several of the study's assumptions are violated by the hardware data (readout errors, crosstalk detected
  on the device alone); the theory reports list them with their consequences.
* **Open questions.** Why the crosstalk detector fires on the device alone; part of the native error that the twin
  does not explain.
* **Statistical reproduction.** Calibrations that sample shot noise reproduce the study's runs in distribution, not draw
  by draw (the study sampled with Aer and its own random streams). Exact results, the hardware analysis and the
  tracking study reproduce to the last digit.
* **Reinforcement learning** is shown from the study's saved results only. The agents are not retrained here, and
  their hyperparameters were not tuned.
