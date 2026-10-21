"""Open Jackson networks of single-server (M/M/1) stations.

Patients arrive at station i from outside at Poisson rate ``external[i]``, are
served at exponential rate ``service[i]``, and then move to station j with
probability ``routing[i][j]`` or leave with probability 1 - sum_j routing[i][j].
Stations are 0-based in code and 1-based in the slides.

``transitions`` and ``inflows`` accept plain numbers or sympy expressions, so
the same code drives the numerical CTMC, the simulation and the symbolic
balance equations in the slides.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


def _shift(state, i, step):
    return tuple(c + step if k == i else c for k, c in enumerate(state))


def transitions(state, external, service, routing):
    """Yield ``(rate, next_state)`` for every transition out of ``state``.

    A component of ``state`` counts as busy unless it is the literal 0, so a
    positive sympy symbol ``n_i`` stands for "at least one patient".
    """
    n = len(state)
    for i in range(n):
        if external[i] != 0:
            yield external[i], _shift(state, i, +1)
    for i in range(n):
        if state[i] == 0:
            continue
        left = _shift(state, i, -1)
        for j in range(n):
            if routing[i][j] != 0:
                yield service[i] * routing[i][j], _shift(left, j, +1)
        leave = 1 - sum(routing[i])
        if abs(leave) > 1e-12:  # a float row such as (1/3, 2/3) leaves round-off
            yield service[i] * leave, left


def _valid(state):
    return all(bool(c >= 0) for c in state)


def inflows(state, external, service, routing):
    """Yield ``(rate, previous_state)`` for every transition into ``state``."""
    n = len(state)
    candidates = []
    for i in range(n):
        candidates.append(_shift(state, i, -1))  # an arrival at i
        candidates.append(_shift(state, i, +1))  # a departure from i
        for j in range(n):
            if j != i:
                candidates.append(_shift(_shift(state, i, +1), j, -1))  # a move i -> j
    seen = set()
    for prev in candidates:
        if prev in seen or not _valid(prev):
            continue
        seen.add(prev)
        for rate, nxt in transitions(prev, external, service, routing):
            if nxt == state:
                yield rate, prev


@dataclass(frozen=True)
class JacksonNetwork:
    external: np.ndarray
    service: np.ndarray
    routing: np.ndarray

    def __post_init__(self):
        external = np.asarray(self.external, dtype=float)
        service = np.asarray(self.service, dtype=float)
        routing = np.asarray(self.routing, dtype=float)
        n = len(service)
        if external.shape != (n,) or routing.shape != (n, n):
            raise ValueError("external, service and routing must describe the same stations")
        if np.any(routing < 0) or np.any(routing.sum(axis=1) > 1 + 1e-9):
            raise ValueError("each routing row must be a sub-probability vector")
        object.__setattr__(self, "external", external)
        object.__setattr__(self, "service", service)
        object.__setattr__(self, "routing", routing)

    @property
    def size(self) -> int:
        return len(self.service)

    def arrival_rates(self) -> np.ndarray:
        """Total arrival rate at each station, from the traffic equations
        lambda_j = lambda_0j + sum_i lambda_i p_ij."""
        return np.linalg.solve(np.eye(self.size) - self.routing.T, self.external)

    def utilisation(self) -> np.ndarray:
        return self.arrival_rates() / self.service

    def is_stable(self) -> bool:
        return bool(np.all(self.utilisation() < 1))

    def _require_stable(self):
        if not self.is_stable():
            rho = np.round(self.utilisation(), 4)
            raise ValueError(f"unstable network: utilisation {rho} has a station at or above 1")

    def mean_number(self) -> np.ndarray:
        """Mean number of patients at each station (waiting or in service)."""
        self._require_stable()
        rho = self.utilisation()
        return rho / (1 - rho)

    def time_per_visit(self) -> np.ndarray:
        """Mean time a patient spends at a station on one visit, 1 / (mu - lambda)."""
        self._require_stable()
        return 1 / (self.service - self.arrival_rates())

    def mean_number_in_system(self) -> float:
        return float(self.mean_number().sum())

    def mean_time_in_system(self) -> float:
        """Mean time from entering to leaving the network, by Little's law."""
        return self.mean_number_in_system() / float(self.external.sum())

    def probability(self, state) -> float:
        """Jackson's product form P(n_1, ..., n_N) = prod_i (1 - rho_i) rho_i^n_i."""
        self._require_stable()
        rho = self.utilisation()
        return float(np.prod((1 - rho) * rho ** np.asarray(state)))
