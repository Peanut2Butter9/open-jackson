"""Numerical stationary distribution of a Jackson network, without the product form.

The network's continuous-time Markov chain is truncated to at most ``cap``
patients per station (transitions that would exceed it are dropped) and its
global balance equations pi Q = 0 are solved directly. As ``cap`` grows the
result converges to Jackson's product form, which makes it an independent check
of both the product form and the balance equations.
"""

from __future__ import annotations

import itertools

import numpy as np
from scipy import sparse

from .network import JacksonNetwork, transitions


def generator(net: JacksonNetwork, cap: int) -> sparse.csr_matrix:
    """Generator matrix Q of the truncated chain, states in row-major order."""
    strides = [(cap + 1) ** (net.size - 1 - k) for k in range(net.size)]
    index = lambda s: sum(c * w for c, w in zip(s, strides))  # noqa: E731 (row-major, like ravel_multi_index)
    external, service, routing = net.external.tolist(), net.service.tolist(), net.routing.tolist()
    rows, cols, rates = [], [], []
    for state in itertools.product(range(cap + 1), repeat=net.size):
        here = index(state)
        for rate, nxt in transitions(state, external, service, routing):
            if max(nxt) <= cap:
                rows.append(here)
                cols.append(index(nxt))
                rates.append(rate)
    size = (cap + 1) ** net.size
    q = sparse.coo_matrix((rates, (rows, cols)), shape=(size, size)).tocsr()
    return q - sparse.diags(np.asarray(q.sum(axis=1)).ravel())


def stationary_distribution(net: JacksonNetwork, cap: int, tol: float = 1e-13, max_steps: int = 10**6) -> np.ndarray:
    """Return pi as an array of shape (cap + 1,) * N, indexed by the state.

    Power iteration on the uniformised chain P = I + Q / q_max, stopped when
    sum |pi Q| < tol. A sparse LU of Q is exact too but fills in badly on this
    3-D lattice (five minutes at cap = 30, against seconds here).
    """
    q = generator(net, cap)
    qt = q.T.tocsr()
    size = q.shape[0]
    step = sparse.identity(size, format="csr") + qt / (-q.diagonal().min() * 1.001)
    pi = np.full(size, 1 / size)
    for k in range(1, max_steps + 1):
        pi = step @ pi
        if k % 500 == 0 and np.abs(qt @ pi).sum() < tol:
            return (pi / pi.sum()).reshape((cap + 1,) * net.size)
    raise RuntimeError(f"power iteration did not converge in {max_steps} steps")


def mean_number(pi: np.ndarray) -> np.ndarray:
    """Mean number of patients at each station under the distribution ``pi``."""
    counts = np.arange(pi.shape[0])
    axes = range(pi.ndim)
    return np.array([(pi.sum(axis=tuple(a for a in axes if a != i)) * counts).sum() for i in axes])
