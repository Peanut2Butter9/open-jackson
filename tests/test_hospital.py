import pytest
import sympy as sp

from openjackson import ctmc, hospital, tex
from openjackson.hospital import lam1, mu1, mu2, mu3, n1, n2, n3


def test_traffic_equations():
    assert hospital.arrival_rates() == [lam1, sp.Rational(13, 15) * lam1, sp.Rational(2, 3) * lam1]


@pytest.mark.parametrize("case", hospital.CASES, ids=str)
def test_product_form_solves_every_balance_equation(case):
    assert hospital.product_form_holds(case)


def test_interior_equation_matches_the_slides():
    out_rate, ins = hospital.balance_equation((n1, n2, n3))
    assert out_rate == lam1 + mu1 + mu2 + mu3
    assert hospital.same_equation(hospital.ORIGINAL_EQUATIONS[(n1, n2, n3)], (out_rate, ins))


@pytest.mark.parametrize("case", [c for c in hospital.CASES if c != (n1, n2, n3)], ids=str)
def test_the_other_original_equations_are_wrong(case):
    original = hospital.ORIGINAL_EQUATIONS[case]
    assert not hospital.same_equation(original, hospital.balance_equation(case))
    assert not hospital.product_form_holds(case, original)


def test_empty_network_balance():
    # Out of (0,0,0) only by an arrival; into it from the doctor and the lab's exits.
    out_rate, ins = hospital.balance_equation((0, 0, 0))
    assert out_rate == lam1
    assert ins == [(mu2, (0, 1, 0)), (mu3 / 5, (0, 0, 1))]


def test_example_is_stable_and_matches_ctmc():
    net = hospital.network()
    assert net.is_stable()
    # Station 1 has the highest load (0.78); 30 patients per station is enough for 2 decimals.
    pi = ctmc.stationary_distribution(net, cap=30)
    assert ctmc.mean_number(pi) == pytest.approx(net.mean_number(), abs=0.02)


def test_generated_latex():
    assert tex.state((n1 - 1, 0, n3 + 1)) == "n_{1}-1,0,n_{3}+1"
    assert tex.rate(sp.Rational(4, 5) * mu3) == r"\tfrac{4}{5}\mu_{3}"
    assert tex.total_rate(mu3 + lam1 + mu2) == r"\lambda_{1} + \mu_{2} + \mu_{3}"
    assert r"\lambda_{1}\,P_{0,0,0} &= \mu_{2}\,P_{0,1,0} + \tfrac{1}{5}\mu_{3}\,P_{0,0,1}" in tex.balance((0, 0, 0))
    assert tex.diagram((n1, n2, n3)).count(r"\draw") == 12
