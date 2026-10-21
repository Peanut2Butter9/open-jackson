import numpy as np
import pytest

from openjackson import paper


@pytest.fixture(scope="module")
def result():
    return paper.reproduce()


def test_table_2_is_counts_over_180_minutes(result):
    assert result["lambda"] == pytest.approx(paper.REPORTED_LAMBDA, abs=0)
    assert result["mu"] == pytest.approx(paper.REPORTED_MU, abs=0)
    # The paper truncates rho: 0.4722 / 0.7722 = 0.61150 is printed as 0.6114.
    assert result["rho"] == pytest.approx(paper.REPORTED_RHO, abs=0)


def test_pseudo_inverse_matches_the_printed_matrix(result):
    # Moore.PNG prints 4 decimals, with some trailing zeros dropped (-0.169).
    assert result["pseudo_inverse_max_abs_diff"] < 3e-4


def test_routing_is_the_normalised_least_squares_solution(result):
    reported = np.array(list(paper.REPORTED_ROUTING.values()))
    assert list(result["routing_normalised"].values()) == pytest.approx(reported, abs=5e-7)
    # Without the normalisation stations 1-3 do not match.
    assert result["routing_least_squares"]["P12"] == pytest.approx(0.4184, abs=1e-4)


def test_routing_system_shape():
    a, b = paper.routing_system(*paper.rates())
    assert a.shape == (11, 9) and b.shape == (11,)
    assert np.linalg.matrix_rank(a) == 8  # the four examination/doctor unknowns are not identifiable


def test_reported_l_is_the_sum_of_rounded_station_values(result):
    assert result["mean_number"] == pytest.approx(37.49, abs=0.01)
    assert result["mean_number_rounded_per_station"] == paper.REPORTED_MEAN_NUMBER


def test_reported_w_adds_up_all_six_stations(result):
    assert result["mean_time_sum_of_stations"] == pytest.approx(paper.REPORTED_MEAN_TIME, abs=1e-3)
    assert result["mean_time_little"] == pytest.approx(36.88, abs=0.01)


def test_reported_routing_overloads_station_5(result):
    rho = result["traffic_utilisation"]
    assert rho[4] > 1
    assert all(r < 1 for k, r in enumerate(rho) if k != 4)
