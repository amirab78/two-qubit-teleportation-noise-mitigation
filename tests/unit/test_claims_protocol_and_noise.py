"""Proven results of steps 1-2 (Parts 9.1-9.2 of the study) as unit tests. Each test names its claim."""
import numpy as np
import pytest

from telecorr import channels as C
from telecorr import protocol as P
from telecorr.paulis import I2
from telecorr.rotations import NOMINAL_NOISE, rotation, u_matrix
from telecorr.states import haar_state, message_state
from telecorr.twirl import message_fidelity

RNG = np.random.default_rng(91)


def random_kraus(rng, n=2):
    """A random one-qubit channel (in general non-unital) with n Kraus operators."""
    G = rng.normal(size=(2 * n, 2)) + 1j * rng.normal(size=(2 * n, 2))
    Q, _ = np.linalg.qr(G)
    return [Q[2 * k:2 * k + 2] for k in range(n)]


def test_T1_5_teleportation_is_the_identity_channel():
    J = P.choi()
    omega = np.eye(4).reshape(16) / 2
    assert np.allclose(J, np.outer(omega, omega), atol=1e-12)
    for _ in range(5):
        psi = haar_state(RNG)
        assert P.fidelity(psi, P.teleport(psi)["bob"]) == pytest.approx(1, abs=1e-12)


def test_T1_6_alice_bits_carry_no_information():
    for _ in range(5):
        A = RNG.normal(size=(4, 4)) + 1j * RNG.normal(size=(4, 4))
        rho = A @ A.conj().T; rho /= np.trace(rho)
        assert np.allclose(P.teleport(rho)["p"], 1 / 16, atol=1e-12)


@pytest.mark.parametrize("frame, value", [("correct", 1.0), ("no correction", 0.25), ("X and Z bits exchanged", 0.4),
                                          ("only X corrections", 0.4), ("only Z corrections", 0.4)])
def test_P1_20_haar_average_of_wrong_wirings(frame, value):
    assert P.average_fidelity(frame=frame) == pytest.approx(value, abs=1e-12)


def test_deferred_and_measured_circuits_are_the_same_channel():
    """The deferred-measurement circuit (Aer, density matrix) equals the exact branch-by-branch simulation."""
    U = C.coherent_unitaries()
    for name in ("product_HT", "bell_30_70"):
        rho_aer = P.aer_bob_state(P.circuit(name, "unitary"), [np.kron(U[3], U[2])])
        rho = P.teleport(message_state(name), legs={2: [U[2]], 3: [U[3]]})["bob"]
        assert np.abs(rho_aer - rho).max() < 1e-12


def test_T2_4_averaging_over_alice_is_a_pauli_twirl():
    """For any channels on q2, q3 (non-unital too) Bob's fidelity is the twirl model's, which needs only diag(R T)."""
    for _ in range(6):
        legs = {2: random_kraus(RNG), 3: random_kraus(RNG, 3)}
        V = {2: u_matrix(*RNG.uniform(-np.pi, np.pi, 3)), 3: u_matrix(*RNG.uniform(-np.pi, np.pi, 3))}
        T = {q: C.T_matrix([V[q] @ K for K in legs[q]]) for q in (2, 3)}
        psi = haar_state(RNG)
        exact = P.fidelity(psi, P.teleport(psi, legs=legs, correction=V)["bob"])
        assert exact == pytest.approx(message_fidelity(psi, {2: np.eye(3), 3: np.eye(3)}, T), abs=1e-12)


def test_T2_5_bob_is_mixed_although_every_branch_is_pure():
    U = C.coherent_unitaries()
    r = P.teleport(message_state("product_HT"), legs={2: [U[2]], 3: [U[3]]})
    assert all(P.purity(b) == pytest.approx(1, abs=1e-12) for b in r["branches"])
    assert P.purity(r["bob"]) < 0.7


def test_P2_6_closed_form_fidelity_of_a_rotation():
    """One leg, residual rotation by alpha about n: F = cos^2(alpha/2) + sin^2(alpha/2) sum_j n_j^2 r_j^2."""
    for _ in range(10):
        n = RNG.normal(size=3); n /= np.linalg.norm(n); alpha = RNG.uniform(0, np.pi)
        a = haar_state(RNG, 2)
        psi = np.kron(np.array([1, 0]), a)                    # message q1 = |0>, q0 = a
        r = [np.real(a.conj() @ s @ a) for s in (C.SIGMA[0], C.SIGMA[1], C.SIGMA[2])]
        legs = {2: [rotation(n, alpha)], 3: [I2]}
        F = P.fidelity(psi, P.teleport(psi, legs=legs)["bob"])
        assert F == pytest.approx(np.cos(alpha / 2) ** 2 + np.sin(alpha / 2) ** 2 * np.dot(n ** 2, np.square(r)), abs=1e-12)


def test_T2_14_no_bob_side_unitary_helps_an_unknown_message():
    U = C.coherent_unitaries()
    legs = {2: [U[2]], 3: [U[3]]}
    base = P.average_fidelity(legs=legs)
    for _ in range(20):
        box = {4: u_matrix(*RNG.uniform(-np.pi, np.pi, 3)), 5: u_matrix(*RNG.uniform(-np.pi, np.pi, 3))}
        assert P.average_fidelity(legs=legs, bob_box=box) <= base + 1e-12


def test_channel_inverse_restores_the_identity():
    U = C.coherent_unitaries()
    V = {q: U[q].conj().T for q in (2, 3)}
    assert P.average_fidelity(legs={2: [U[2]], 3: [U[3]]}, correction=V) == pytest.approx(1, abs=1e-12)


def test_noise_model_is_valid_and_matches_its_closed_form():
    for fam in C.FAMILIES:
        legs = C.leg_kraus(None, fam, 0.2)
        T = C.leg_T(None, fam, 0.2)
        for q in (2, 3):
            assert C.is_cptp(legs[q])
            assert np.allclose(C.T_matrix(legs[q]), T[q], atol=1e-12)
    assert not C.is_unital(C.incoherent_kraus(amp=0.3)) and C.is_unital(C.incoherent_kraus(depol=0.3))
    assert NOMINAL_NOISE[2]["angle"] == 0.5 and NOMINAL_NOISE[3]["angle"] == 1.3
