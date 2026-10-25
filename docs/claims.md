# Proven results and the tests that check them

The research project labels its results **T** (theorem: proven in general), **P** (proposition: proven under stated
assumptions) and **O** (observation: data or simulation); the numbering follows the study's theory reports (Part 9.k
for step k). Every proven result used by this repository is a unit test, and the saved numbers of the study are
regression tests. Run all of them with `pytest` (about 1.5 minutes; `pytest -m "not slow"` skips the one-minute
tracking test).

## Proven results (unit tests)

| claim | what is checked | test |
|---|---|---|
| T1.5 | Teleportation is the identity channel | `tests/unit/test_claims_protocol_and_noise.py::test_T1_5_teleportation_is_the_identity_channel` |
| T1.6 | Alice bits carry no information | `tests/unit/test_claims_protocol_and_noise.py::test_T1_6_alice_bits_carry_no_information` |
| P1.20 | Haar average of wrong wirings | `tests/unit/test_claims_protocol_and_noise.py::test_P1_20_haar_average_of_wrong_wirings` |
| T2.4 | For any channels on q2, q3 (non-unital too) Bob's fidelity is the twirl model's, which needs only diag(R T). | `tests/unit/test_claims_protocol_and_noise.py::test_T2_4_averaging_over_alice_is_a_pauli_twirl` |
| T2.5 | Bob is mixed although every branch is pure | `tests/unit/test_claims_protocol_and_noise.py::test_T2_5_bob_is_mixed_although_every_branch_is_pure` |
| P2.6 | One leg, residual rotation by alpha about n: F = cos^2(alpha/2) + sin^2(alpha/2) sum_j n_j^2 r_j^2. | `tests/unit/test_claims_protocol_and_noise.py::test_P2_6_closed_form_fidelity_of_a_rotation` |
| T2.14 | No bob side unitary helps an unknown message | `tests/unit/test_claims_protocol_and_noise.py::test_T2_14_no_bob_side_unitary_helps_an_unknown_message` |
| T3.5 | All four outcome probabilities of every probe circuit: twirl formula = the qiskit circuit (Aer, density matrix). | `tests/unit/test_claims_correction.py::test_T3_5_probe_probabilities_are_products_of_leg_factors` |
| T3.6, T5.9 | The cost is sandwiched by the average fidelity | `tests/unit/test_claims_correction.py::test_T3_6_T5_9_the_cost_is_sandwiched_by_the_average_fidelity` |
| T3.8 | One probe is blind to rotations about its own axis | `tests/unit/test_claims_correction.py::test_T3_8_one_probe_is_blind_to_rotations_about_its_own_axis` |
| P4.1 | Cost is bilinear in the squared quaternions | `tests/unit/test_claims_correction.py::test_P4_1_cost_is_bilinear_in_the_squared_quaternions` |
| T4.2 | Along the path u(t) = (1 - t) u + t e0 (both legs) the cost decreases strictly to 0 from every point. | `tests/unit/test_claims_correction.py::test_T4_2_no_spurious_local_minimum` |
| T4.3 | Critical values | `tests/unit/test_claims_correction.py::test_T4_3_critical_values` |
| P4.10 | An exact cost C < 0.01 implies F_e > 0.985, a lower bound on the fidelity of every message. | `tests/unit/test_claims_correction.py::test_P4_10_an_accepted_start_guarantees_a_high_fidelity` |
| T5.4 | Average fidelity formula | `tests/unit/test_claims_correction.py::test_T5_4_average_fidelity_formula` |
| T5.6 | Ceiling is procrustes and the inverse for symmetric noise | `tests/unit/test_claims_correction.py::test_T5_6_ceiling_is_procrustes_and_the_inverse_for_symmetric_noise` |
| P5.10 | Amplitude damping p = 0.3: the 3-probe optimum misses the ceiling by about 1e-6 (values of the study). | `tests/unit/test_claims_correction.py::test_P5_10_probe_bias_is_real_but_second_order` |
| P5.15 | Rotosolve coordinates are exact sinusoids | `tests/unit/test_claims_correction.py::test_P5_15_rotosolve_coordinates_are_exact_sinusoids` |
| T6.2, T6.3 | Two qubit twirl and fidelities | `tests/unit/test_claims_correlated_and_tracking.py::test_T6_2_T6_3_two_qubit_twirl_and_fidelities` |
| T6.5 | Inverse is optimal for symmetric psd noise | `tests/unit/test_claims_correlated_and_tracking.py::test_T6_5_inverse_is_optimal_for_symmetric_psd_noise` |
| T6.7 | Pure zz crosstalk local ceiling | `tests/unit/test_claims_correlated_and_tracking.py::test_T6_7_pure_zz_crosstalk_local_ceiling` |
| T6.11 | Tomo75 reads every entry exactly once | `tests/unit/test_claims_correlated_and_tracking.py::test_T6_11_tomo75_reads_every_entry_exactly_once` |
| P6.13 | Detector signal is zero for product channels | `tests/unit/test_claims_correlated_and_tracking.py::test_P6_13_detector_signal_is_zero_for_product_channels` |
| T7.2 | Tomo9 reads every entry of t exactly once | `tests/unit/test_claims_correlated_and_tracking.py::test_T7_2_tomo9_reads_every_entry_of_T_exactly_once` |
| T7.3 | Procrustes batch equals the closed form value | `tests/unit/test_claims_correlated_and_tracking.py::test_T7_3_procrustes_batch_equals_the_closed_form_value` |
| P7.8 | Sequential kalman updates equal the batch information form | `tests/unit/test_claims_correlated_and_tracking.py::test_P7_8_sequential_kalman_updates_equal_the_batch_information_form` |
| P7.9 | Steady state of a random walk filter | `tests/unit/test_claims_correlated_and_tracking.py::test_P7_9_steady_state_of_a_random_walk_filter` |
| T8.1 | T8.1: for every basis B and correction Z^z X^x, (Z^z X^x)^dag Pi_b^B (Z^z X^x) = Pi_{b xor f}^B with | `tests/unit/test_hardware_claims.py::test_T8_1_pauli_frame_equals_live_correction` |
| T8.1 | Alice's bit 0 / 1 = Z corrections of q4 / q5, bits 2 / 3 = X corrections (Bob key = b5 * 2 + b4). | `tests/unit/test_hardware_claims.py::test_T8_1_frame_keys_bit_layout` |
| T8.2 | T8.2: a bit with flip probabilities e0 = P(1|0), e1 = P(0|1) measures lambda_hat = (1 - e0 - e1) lambda + (e1 - e0). | `tests/unit/test_hardware_claims.py::test_T8_2_readout_law` |
| T8.3 | T8.3: with exact probabilities the linear estimate returns the state, so F_lin is unbiased. | `tests/unit/test_hardware_claims.py::test_T8_3_linear_tomography_is_exact_on_probabilities` |
| T8.3 | F_lin = sum w f (the weights reproduce tomography's F on any counts); its multinomial variance formula agrees | `tests/unit/test_hardware_claims.py::test_T8_3_linear_F_is_linear_in_frequencies_and_variance_exact` |
| T8.3 | Mle is a density matrix and fixes physical states | `tests/unit/test_hardware_claims.py::test_T8_3_mle_is_a_density_matrix_and_fixes_physical_states` |
| P8.10 | P8.10: Var f = (s_c^2 + f^2 s_a^2 + (1 - f)^2 s_n^2) / D^2 for f = (c - n)/(a - n) (numeric derivative check), | `tests/unit/test_hardware_claims.py::test_P8_10_delta_method_and_fieller` |

## Further unit tests

| what is checked | test |
|---|---|
| Rotosolve with exact costs descends monotonically | `tests/unit/test_claims_correction.py::test_rotosolve_with_exact_costs_descends_monotonically` |
| Recentre keeps the correction unchanged at zero | `tests/unit/test_claims_correction.py::test_recentre_keeps_the_correction_unchanged_at_zero` |
| Zz is the pauli used for correlated dephasing | `tests/unit/test_claims_correlated_and_tracking.py::test_zz_is_the_pauli_used_for_correlated_dephasing` |
| The deferred-measurement circuit (Aer, density matrix) equals the exact branch-by-branch simulation. | `tests/unit/test_claims_protocol_and_noise.py::test_deferred_and_measured_circuits_are_the_same_channel` |
| Channel inverse restores the identity | `tests/unit/test_claims_protocol_and_noise.py::test_channel_inverse_restores_the_identity` |
| Noise model is valid and matches its closed form | `tests/unit/test_claims_protocol_and_noise.py::test_noise_model_is_valid_and_matches_its_closed_form` |
| confusion_from_counts recovers e0, e1 and Bob's matrix from pattern counts built with known errors. | `tests/unit/test_hardware_claims.py::test_readout_calibration_counts` |

## Saved results of the study (regression tests)

| what is reproduced | test |
|---|---|
| Manifest hashes | `tests/regression/test_data_integrity.py::test_manifest_hashes` |
| Predictions hashed before sending | `tests/regression/test_data_integrity.py::test_predictions_hashed_before_sending` |
| approved hash = frozen job hash = hash recorded when the job was sent, for all 4 jobs. | `tests/regression/test_data_integrity.py::test_job_hash_chain` |
| Counts are complete | `tests/regression/test_data_integrity.py::test_counts_are_complete` |
| No account information in data | `tests/regression/test_data_integrity.py::test_no_account_information_in_data` |
| Every prediction file of the replication has the SHA-256 recorded before sending, the frozen job names that | `tests/regression/test_data_integrity.py::test_replication_predictions_hashed_before_sending` |
| approved hash = frozen job hash = hash recorded when the job was sent, and usage <= the approved cap. | `tests/regression/test_data_integrity.py::test_replication_job_hash_chain` |
| Replication counts are complete | `tests/regression/test_data_integrity.py::test_replication_counts_are_complete` |
| Round1 equals saved | `tests/regression/test_hardware.py::test_round1_equals_saved` |
| Round2 equals saved | `tests/regression/test_hardware.py::test_round2_equals_saved` |
| Corrections and n0 rule | `tests/regression/test_hardware.py::test_corrections_and_n0_rule` |
| Per circuit checks | `tests/regression/test_hardware.py::test_per_circuit_checks` |
| Compare and hypotheses | `tests/regression/test_hardware.py::test_compare_and_hypotheses` |
| Twirl model equals preregistered file | `tests/regression/test_hardware.py::test_twirl_model_equals_preregistered_file` |
| All saved results match | `tests/regression/test_hardware.py::test_all_saved_results_match` |
| The numbers quoted in the paper / README (rounded as printed there). | `tests/regression/test_hardware.py::test_headline_numbers` |
| All saved results match | `tests/regression/test_replication.py::test_all_saved_results_match` |
| Round1 equals saved | `tests/regression/test_replication.py::test_round1_equals_saved` |
| K2 / F2 ran the corrections that the prediction files recorded before sending, and those are the | `tests/regression/test_replication.py::test_round2_was_built_from_round1` |
| Three weeks later on the same device: the fresh correction recovers >= 90 % of N1's loss and reaches the | `tests/regression/test_replication.py::test_replication_on_kingston` |
| ibm_fez, new chain, same pipeline: N1 0.706 / 0.511 -> 0.910 / 0.891; both pre-registered criteria pass. | `tests/regression/test_replication.py::test_second_device_without_tuning` |
| The waiting phase grows linearly (axis -z, both legs near 15 kHz today); a fresh TOMO-9 recovers more than | `tests/regression/test_replication.py::test_waits_linear_and_corrected` |
| Part 9.8's spectator phase holds as a DIFFERENCE (K3 - K4 = pi zeta t) but not absolutely: without spectator | `tests/regression/test_replication.py::test_dd_spectator_difference_and_absolute_failure` |
| In one job on Fez: product_HT's 'corrected above native' was drift, bell_30_70's is real, and removing Fez's | `tests/regression/test_replication.py::test_g1_same_job` |
| z(t) = a + b t under DD: the slope is about 8 sigma from zero on both legs (+2.5 / +2.7 kHz), so CONST is | `tests/regression/test_replication.py::test_g2_dd_rotation_grows_with_the_wait` |
| G3-A fails at both waits; DD + correction is worse than K2's fresh correction without DD by > 4 sigma. | `tests/regression/test_replication.py::test_g3_correction_under_dd_barely_helps` |
| Jobs within their caps | `tests/regression/test_replication.py::test_jobs_within_their_caps` |
| Step1 unpreparation reads 00 in every shot | `tests/regression/test_simulation.py::test_step1_unpreparation_reads_00_in_every_shot` |
| Step2 fidelity and purity of the coherent channel | `tests/regression/test_simulation.py::test_step2_fidelity_and_purity_of_the_coherent_channel` |
| Step2 best bob side correction reaches lambda max | `tests/regression/test_simulation.py::test_step2_best_bob_side_correction_reaches_lambda_max` |
| Stage 3's saved angles: probe cost 8.93e-5, F = 0.999906 / 0.999874, residual 0.0164 / 0.0163 rad. | `tests/regression/test_simulation.py::test_step3_saved_stage3_correction` |
| A fresh COBYLA calibration from counts (seed fixed): F_avg > 0.995 (30 seeds: median 1 - F_avg = 4.6e-4). | `tests/regression/test_simulation.py::test_step3_fresh_calibration_reaches_the_shot_noise_floor` |
| Step4 hard case | `tests/regression/test_simulation.py::test_step4_hard_case` |
| Step5 realistic channel e4 | `tests/regression/test_simulation.py::test_step5_realistic_channel_E4` |
| Step6 nominal correlated channel table | `tests/regression/test_simulation.py::test_step6_nominal_correlated_channel_table` |
| Stage 7f, scenario 'slow' (200 episodes, its seeds): ABC and TOMO-9 every step reproduce the saved values. | `tests/regression/test_simulation.py::test_step7_tracking_slow_drift_exactly` |
