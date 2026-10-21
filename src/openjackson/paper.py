"""Reproduce the six-station hospital network of Putri and Endrayanto,
"On the Application of the Open Jackson Queuing Network in Hospital".

Stations: 1 registration, 2 and 3 examination rooms, 4 and 5 doctor's rooms,
6 administration and pharmacy. Patients go 1 -> {2, 3} -> {4, 5} -> 6 -> out.
The ``REPORTED_*`` values are copied from the paper's tables and conclusion
as quoted in the slides (Tabel1.png, Tabel2.png, Moore.PNG, main.tex).
"""

from __future__ import annotations

import numpy as np

from .network import JacksonNetwork

ARRIVALS = np.array([183, 85, 110, 85, 78, 195])
DEPARTURES = np.array([235, 139, 146, 130, 86, 206])
# Not stated in the slides; every lambda and mu in the paper's Table 2 is the
# count divided by 180, i.e. a three-hour observation in minutes.
WINDOW_MINUTES = 180

REPORTED_LAMBDA = np.array([1.0167, 0.4722, 0.6111, 0.4722, 0.4333, 1.0833])
REPORTED_MU = np.array([1.3056, 0.7722, 0.8111, 0.7222, 0.4778, 1.1444])
REPORTED_RHO = np.array([0.7787, 0.6114, 0.7534, 0.6538, 0.9068, 0.9466])
REPORTED_ROUTING = {
    "P12": 0.4436053, "P13": 0.5563947, "P24": 0.5153421, "P25": 0.4846579, "P34": 0.5159005,
    "P35": 0.4840995, "P46": 0.9461371, "P56": 0.9185853, "P60": 1.0000004,
}  # fmt: skip
REPORTED_MEAN_NUMBER = 39  # patients in the whole network
REPORTED_MEAN_TIME = 54.6101  # minutes in the whole network

# Unknowns of the routing system, in the paper's order; station 0 is "outside".
EDGES = [(1, 2), (1, 3), (2, 4), (2, 5), (3, 4), (3, 5), (4, 6), (5, 6), (6, 0)]
EDGE_NAMES = [f"P{i}{j}" for i, j in EDGES]

# The pseudo-inverse printed in the slides (Moore.PNG): 9 x 11, columns in the
# order of the right-hand side (lambda_2..lambda_6, mu_1..mu_6).
REPORTED_PSEUDO_INVERSE = np.array([
    [0.5106, -0.2553, 0, 0, 0, 0.2553, 0, 0, 0, 0, 0],
    [-0.2553, 0.5106, 0, 0, 0, 0.2553, 0, 0, 0, 0, 0],
    [0, 0, 0.4697, -0.1459, 0, 0, 0.4856, -0.1619, 0, 0, 0],
    [0, 0, -0.1459, 0.4697, 0, 0, 0.4856, -0.1619, 0, 0, 0],
    [0, 0, 0.4775, -0.169, 0, 0, -0.1541, 0.4623, 0, 0, 0],
    [0, 0, -0.1692, 0.4775, 0, 0, -0.1541, 0.4623, 0, 0, 0],
    [0, 0, 0, 0, 0.4615, 0, 0, 0, 0.9231, -0.4615, 0],
    [0, 0, 0, 0, 0.6976, 0, 0, 0, -0.6976, 1.3953, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.8738],
])  # fmt: skip


def rates() -> tuple[np.ndarray, np.ndarray]:
    """lambda_i and mu_i in patients per minute, rounded to 4 decimals as in the paper."""
    return np.round(ARRIVALS / WINDOW_MINUTES, 4), np.round(DEPARTURES / WINDOW_MINUTES, 4)


def utilisation(lam, mu) -> np.ndarray:
    """rho_i = lambda_i / mu_i, truncated (not rounded) to 4 decimals as the paper's table is."""
    return np.floor(lam / mu * 1e4) / 1e4


def routing_system(lam, mu) -> tuple[np.ndarray, np.ndarray]:
    """The paper's 11 x 9 linear system A p = b for the routing probabilities.

    Rows 1-5: lambda_j = sum_i mu_i P_ij for j = 2..6 (flow into each station).
    Rows 6-11: mu_i = sum_j mu_i P_ij for i = 1..6 (each station's outflow is split).
    The slide's copy of the system (Balance Eqs.png) also puts mu_3 P_35 in the
    lambda_6 row, but its printed pseudo-inverse was computed without it.
    """
    a = np.zeros((11, len(EDGES)))
    b = np.concatenate([lam[1:], mu])
    for k, (i, j) in enumerate(EDGES):
        if j:
            a[j - 2, k] = mu[i - 1]
        a[5 + i - 1, k] = mu[i - 1]
    return a, b


def routing_least_squares(lam, mu) -> np.ndarray:
    """Minimum-norm least-squares solution pinv(A) b, in ``EDGES`` order."""
    a, b = routing_system(lam, mu)
    return np.linalg.pinv(a) @ b


def normalise_branches(p) -> np.ndarray:
    """Rescale the two branches out of stations 1, 2 and 3 to sum to 1.

    This step is not in the slides, but it turns the least-squares solution into
    the reported routing to 7 decimals. Stations 4, 5 and 6 have a single
    successor and keep their least-squares value.
    """
    p = np.asarray(p, dtype=float).copy()
    for station in (1, 2, 3):
        idx = [k for k, (i, _) in enumerate(EDGES) if i == station]
        p[idx] /= p[idx].sum()
    return p


def mean_number_per_station(rho) -> np.ndarray:
    return rho / (1 - rho)


def network_from_routing(routing: dict, lam1: float, mu) -> JacksonNetwork:
    """The Jackson network implied by external arrivals at station 1 and ``routing``."""
    p = np.zeros((6, 6))
    for name, value in routing.items():
        i, j = int(name[1]), int(name[2])
        if j:
            p[i - 1, j - 1] = value
    return JacksonNetwork(external=[lam1, 0, 0, 0, 0, 0], service=mu, routing=p)


def reproduce() -> dict:
    """Every number of the paper next to its reproduction."""
    lam, mu = rates()
    rho = utilisation(lam, mu)
    a, _ = routing_system(lam, mu)
    least_squares = routing_least_squares(lam, mu)
    normalised = normalise_branches(least_squares)
    per_station = mean_number_per_station(rho)
    net = network_from_routing(dict(zip(EDGE_NAMES, normalised)), lam[0], mu)
    return {
        "lambda": lam.tolist(),
        "mu": mu.tolist(),
        "rho": rho.tolist(),
        "rho_exact": (ARRIVALS / DEPARTURES).tolist(),
        "pseudo_inverse_max_abs_diff": float(np.abs(np.linalg.pinv(a) - REPORTED_PSEUDO_INVERSE).max()),
        "routing_least_squares": dict(zip(EDGE_NAMES, least_squares.tolist())),
        "routing_normalised": dict(zip(EDGE_NAMES, normalised.tolist())),
        "routing_max_abs_diff": float(np.abs(normalised - np.array(list(REPORTED_ROUTING.values()))).max()),
        "mean_number_per_station": per_station.tolist(),
        "mean_number": float(per_station.sum()),
        "mean_number_rounded_per_station": int(np.round(per_station).sum()),
        # The paper's waiting time: the per-visit times of all six stations added
        # up, as if every patient visited every station.
        "mean_time_sum_of_stations": float((per_station / lam).sum()),
        # Little's law: mean number in the network over the rate patients enter it.
        "mean_time_little": float(per_station.sum() / lam[0]),
        "traffic_arrival_rates": net.arrival_rates().tolist(),
        "traffic_utilisation": net.utilisation().tolist(),
    }
