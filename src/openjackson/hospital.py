"""The slides' three-station model: registration (1), doctor (2), laboratory (3).

Patients arrive at registration at rate lambda_1. After registration a third go
to the doctor and two thirds to the laboratory. After the laboratory four
fifths go on to the doctor and one fifth leave. Everyone leaves after the doctor.
"""

from __future__ import annotations

import sympy as sp

from .network import JacksonNetwork, inflows, transitions

lam1, mu1, mu2, mu3 = sp.symbols("lambda_1 mu_1 mu_2 mu_3", positive=True)
n1, n2, n3 = sp.symbols("n_1 n_2 n_3", positive=True, integer=True)
RATES = (lam1, mu1, mu2, mu3)

EXTERNAL = (lam1, 0, 0)
SERVICE = (mu1, mu2, mu3)
ROUTING = (
    (0, sp.Rational(1, 3), sp.Rational(2, 3)),
    (0, 0, 0),
    (0, sp.Rational(4, 5), 0),
)

# The eight balance-equation cases, in the order of the slides: a symbol n_i
# means station i holds at least one patient, 0 means it is empty.
CASES = [
    (0, 0, 0),
    (n1, 0, 0),
    (0, n2, 0),
    (0, 0, n3),
    (n1, n2, n3),
    (n1, n2, 0),
    (0, n2, n3),
    (n1, 0, n3),
]

# The balance equations as the original slides write them, (rate out, [(rate,
# previous state), ...]) per case. Case 2 has no "=" and is read as
# rate-out-term = the other terms; a bare 2/3, 1/5 or mu is read as the
# nearest rate (2/3 mu_1, 1/5 mu_3, mu_3) to give the slides the benefit of the doubt.
_t, _f = sp.Rational(1, 3), sp.Rational(1, 5)
ORIGINAL_EQUATIONS = {
    (0, 0, 0): (lam1, [(_t * mu1, (0, 1, 0)), (2 * _t * mu1, (0, 0, 1))]),
    (n1, 0, 0): (lam1 + mu1, [(_f * mu3, (n1, 0, 1)), (mu2, (n1, 1, 0))]),
    (0, n2, 0): (lam1 + mu1, [(_t * mu1, (1, n2 - 1, 0)), (mu2, (0, n2 + 1, 0)), (4 * _f * mu3, (0, n2 - 1, 1))]),
    (0, 0, n3): (lam1 + mu3, [(2 * _t * mu1, (1, 0, n3 - 1)), (mu3, (0, 0, n3 + 1))]),
    (n1, n2, n3): (
        lam1 + mu1 + mu2 + mu3,
        [
            (lam1, (n1 - 1, n2, n3)),
            (_t * mu1, (n1 + 1, n2 - 1, n3)),
            (2 * _t * mu1, (n1 + 1, n2, n3 - 1)),
            (mu2, (n1, n2 + 1, n3)),
            (_f * mu3, (n1, n2, n3 + 1)),
            (4 * _f * mu3, (n1, n2 - 1, n3 + 1)),
        ],
    ),
    (n1, n2, 0): (2 * _t * mu1 + lam1 + mu2, [(_t * mu1, (n1 + 1, n2 - 1, 0)), (4 * _f * mu3, (n1, n2 - 1, 1))]),
    (0, n2, n3): (
        mu3 + mu2 + lam1,
        [
            (2 * _t * mu1, (1, n2, n3 + 1)),
            (_f * mu3, (0, n2, n3 + 1)),
            (mu2, (0, n2 + 1, n3)),
            (_t * mu1, (1, n2 + 1, n3)),
        ],
    ),
    (n1, 0, n3): (_t * mu1 + mu3 + lam1, [(2 * _t * mu1, (n1 - 1, 0, n3 + 1))]),
}

# Illustrative rates in patients per minute (the slides give none). Registration
# is the paper's station 1; the doctor pools its two doctor's rooms (stations 4
# and 5) and the laboratory its two examination rooms (stations 2 and 3).
EXAMPLE = {lam1: 1.0167, mu1: 1.3056, mu2: 0.7222 + 0.4778, mu3: 0.7722 + 0.8111}


def rate_key(rate):
    """Sort key that lists rates as lambda_1, mu_1, mu_2, mu_3, smaller share first."""
    coef, symbol = rate.as_coeff_Mul()
    return RATES.index(symbol), coef


def out_transitions(state):
    return sorted(transitions(state, EXTERNAL, SERVICE, ROUTING), key=lambda t: rate_key(t[0]))


def in_transitions(state):
    return sorted(inflows(state, EXTERNAL, SERVICE, ROUTING), key=lambda t: rate_key(t[0]))


def balance_equation(state):
    """(total rate out of ``state``, [(rate, previous state), ...]) of its global balance equation."""
    return sum(rate for rate, _ in out_transitions(state)), in_transitions(state)


def arrival_rates():
    """Symbolic traffic equations: Lambda = lambda_0 + P^T Lambda."""
    p = sp.Matrix(ROUTING)
    return list((sp.eye(3) - p.T).solve(sp.Matrix(EXTERNAL)))


def utilisation():
    return [lam / mu for lam, mu in zip(arrival_rates(), SERVICE)]


def same_equation(a, b) -> bool:
    """Do two balance equations (rate out, [(rate, previous state), ...]) agree term by term?"""
    return sp.simplify(a[0] - b[0]) == 0 and sorted(map(str, a[1])) == sorted(map(str, b[1]))


def product_form_holds(state, equation=None) -> bool:
    """Does P(n) = C prod rho_i^n_i satisfy the balance equation of ``state``?

    Dividing both sides by P(state) leaves prod rho_i^(m_i - n_i) for each
    neighbour m, so the check is an identity in lambda_1 and the mu_i.
    ``equation`` defaults to the one derived from the transitions.
    """
    rho = utilisation()
    out_rate, ins = equation or balance_equation(state)
    inflow = sum(
        rate * sp.Mul(*(r ** sp.simplify(m - s) for r, m, s in zip(rho, prev, state))) for rate, prev in ins
    )
    return sp.simplify(out_rate - inflow) == 0


def network(rates=EXAMPLE) -> JacksonNetwork:
    return JacksonNetwork(
        external=[float(sp.sympify(e).subs(rates)) for e in EXTERNAL],
        service=[float(s.subs(rates)) for s in SERVICE],
        routing=[[float(p) for p in row] for row in ROUTING],
    )
