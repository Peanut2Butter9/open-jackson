import numpy as np
import pytest

from openjackson import JacksonNetwork, ctmc, inflows, transitions
from openjackson.simulate import simulate


def tandem(lam=1.0, mu1=2.0, mu2=3.0):
    return JacksonNetwork(external=[lam, 0], service=[mu1, mu2], routing=[[0, 1], [0, 0]])


def test_single_mm1_matches_textbook():
    net = JacksonNetwork(external=[1.0], service=[4.0], routing=[[0.0]])
    assert net.utilisation() == pytest.approx([0.25])
    assert net.mean_number() == pytest.approx([1 / 3])  # rho / (1 - rho)
    assert net.time_per_visit() == pytest.approx([1 / 3])  # 1 / (mu - lambda)
    assert net.probability([2]) == pytest.approx(0.75 * 0.25**2)


def test_traffic_equations_with_feedback():
    # Station 1 sends half its patients back to itself: lambda_1 = 1 + lambda_1 / 2.
    net = JacksonNetwork(external=[1.0], service=[5.0], routing=[[0.5]])
    assert net.arrival_rates() == pytest.approx([2.0])


def test_little_law_in_a_tandem():
    net = tandem()
    assert net.mean_time_in_system() == pytest.approx(net.time_per_visit().sum())


def test_unstable_network_is_rejected():
    net = tandem(lam=2.5)
    assert not net.is_stable()
    with pytest.raises(ValueError, match="unstable"):
        net.mean_number()


def test_routing_rows_must_be_sub_probabilities():
    with pytest.raises(ValueError):
        JacksonNetwork(external=[1, 0], service=[2, 2], routing=[[0, 1.2], [0, 0]])


def test_inflows_invert_transitions():
    net = tandem()
    args = (net.external, net.service, net.routing)
    state = (2, 1)
    for rate, prev in inflows(state, *args):
        assert (rate, state) in list(transitions(prev, *args))
    # every state that can reach (2, 1) in one step: (1, 1) by an arrival,
    # (3, 0) by a move 1 -> 2, and (2, 2) by a departure from station 2
    assert sorted(prev for _, prev in inflows(state, *args)) == [(1, 1), (2, 2), (3, 0)]


def test_truncated_ctmc_converges_to_product_form():
    net = tandem(lam=1.0, mu1=3.0, mu2=4.0)
    pi = ctmc.stationary_distribution(net, cap=25)
    assert pi.sum() == pytest.approx(1)
    assert ctmc.mean_number(pi) == pytest.approx(net.mean_number(), rel=1e-6)
    assert pi[1, 2] == pytest.approx(net.probability((1, 2)), rel=1e-6)


def test_simulation_agrees_with_formulas():
    net = tandem(lam=1.0, mu1=2.0, mu2=2.5)
    result = simulate(net, horizon=40_000, seed=1)
    assert result.mean_number == pytest.approx(net.mean_number(), rel=0.1)
    assert result.mean_time_in_system == pytest.approx(net.mean_time_in_system(), rel=0.1)
    assert np.all(result.mean_number > 0)
