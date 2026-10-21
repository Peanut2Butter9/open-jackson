"""Patient-level discrete-event simulation of a Jackson network.

Each station serves one patient at a time, first come first served. The
simulation measures the time-average number of patients at each station and
the mean time each patient spends in the network. It shares no code with the
analytical formulas, so it checks them, Little's law included.
"""

from __future__ import annotations

import bisect
import heapq
from collections import deque
from dataclasses import dataclass
from itertools import accumulate

import numpy as np

from .network import JacksonNetwork


@dataclass(frozen=True)
class SimulationResult:
    mean_number: np.ndarray  # time-average number of patients at each station
    mean_time_in_system: float  # mean time from entering to leaving the network
    patients: int  # patients that entered after the warm-up and have left


def simulate(net: JacksonNetwork, horizon: float, seed: int = 0, warmup: float = 0.1) -> SimulationResult:
    """Simulate for ``horizon`` time units and discard the first ``warmup`` fraction."""
    rng = np.random.default_rng(seed)
    n = net.size
    start = warmup * horizon
    mean_service = (1 / net.service).tolist()
    # Cumulative routing per station; drawing past the last entry means leaving.
    cum = [list(accumulate(row)) for row in net.routing.tolist()]
    queues = [deque() for _ in range(n)]  # entry times of the patients at each station
    events = []  # (time, kind, station); kind 0 = external arrival, 1 = service completion
    for i in range(n):
        if net.external[i] > 0:
            heapq.heappush(events, (rng.exponential(1 / net.external[i]), 0, i))

    area = [0.0] * n
    last = start
    sojourn_sum, patients = 0.0, 0

    def join(i, entered, now):
        queues[i].append(entered)
        if len(queues[i]) == 1:
            heapq.heappush(events, (now + rng.exponential(mean_service[i]), 1, i))

    while events:
        now, kind, i = heapq.heappop(events)
        if now > horizon:
            break
        if now > start:
            for k in range(n):
                area[k] += len(queues[k]) * (now - last)
            last = now
        if kind == 0:
            join(i, now, now)
            heapq.heappush(events, (now + rng.exponential(1 / net.external[i]), 0, i))
            continue
        entered = queues[i].popleft()
        if queues[i]:
            heapq.heappush(events, (now + rng.exponential(mean_service[i]), 1, i))
        j = bisect.bisect_right(cum[i], rng.random())
        if j < n:
            join(j, entered, now)
        elif entered > start:
            sojourn_sum += now - entered
            patients += 1

    for k in range(n):
        area[k] += len(queues[k]) * (horizon - last)
    return SimulationResult(np.array(area) / (horizon - start), sojourn_sum / max(patients, 1), patients)
