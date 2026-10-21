"""LaTeX snippets for the slides, generated from the model so they cannot drift from it."""

from __future__ import annotations

import sympy as sp

from . import hospital


def rate(r) -> str:
    coef, symbol = r.as_coeff_Mul()
    head = "" if coef == 1 else rf"\tfrac{{{coef.p}}}{{{coef.q}}}"
    return head + sp.latex(symbol)


def total_rate(r) -> str:
    terms = sorted(sp.Add.make_args(sp.expand(r)), key=hospital.rate_key)
    return " + ".join(rate(t) for t in terms)


def state(s) -> str:
    return ",".join(str(c) if isinstance(c, int) else sp.latex(c).replace(" ", "") for c in s)


def prob(s) -> str:
    return f"P_{{{state(s)}}}"


def condition(s) -> str:
    """The cases' condition, e.g. n_1>0, n_2=0, n_3>0."""
    return ",\\ ".join(f"n_{{{i}}}{'=0' if c == 0 else '>0'}" for i, c in enumerate(s, start=1))


def balance(s, per_line: int = 3) -> str:
    """The global balance equation of case ``s`` as a display."""
    out_rate, ins = hospital.balance_equation(s)
    total = total_rate(out_rate)
    lhs = f"({total})\\,{prob(s)}" if "+" in total else f"{total}\\,{prob(s)}"
    terms = [f"{rate(r)}\\,{prob(prev)}" for r, prev in ins]
    lines = [" + ".join(terms[k : k + per_line]) for k in range(0, len(terms), per_line)]
    body = f"{lhs} &= {lines[0]}" + "".join(f"\\\\\n  &\\quad + {line}" for line in lines[1:])
    return f"\\[\n\\begin{{aligned}}\n  {body}\n\\end{{aligned}}\n\\]\n"


def transition_rows() -> str:
    """Rows of the transition table: from, to, rate, condition."""
    generic = (hospital.n1, hospital.n2, hospital.n3)
    rows = []
    for r, nxt in hospital.out_transitions(generic):
        coef, symbol = r.as_coeff_Mul()
        need = "--" if symbol == hospital.lam1 else f"$n_{{{hospital.SERVICE.index(symbol) + 1}}}>0$"
        rows.append(f"$({state(generic)})$ & $({state(nxt)})$ & ${rate(r)}$ & {need} \\\\")
    return "\n".join(rows) + "\n"


def diagram(s, dx: float = 4.6, dy: float = 0.95) -> str:
    """TikZ state-transition diagram of case ``s``: inflows on the left, outflows on the right."""
    ins, outs = hospital.in_transitions(s), hospital.out_transitions(s)
    lines = [
        r"\begin{tikzpicture}[>=stealth, font=\small,",
        r"  st/.style={draw, rounded corners=7pt, minimum height=6mm, inner sep=3pt},",
        r"  into/.style={st, draw=mygreen, text=mygreen},",
        r"  rt/.style={fill=white, inner sep=1pt, font=\footnotesize}]",
        rf"\node[st, very thick, font=\normalsize] (c) at (0,0) {{${state(s)}$}};",
    ]
    for k, (r, prev) in enumerate(ins):
        y = ((len(ins) - 1) / 2 - k) * dy
        lines.append(rf"\node[into] (i{k}) at ({-dx},{y:.3f}) {{${state(prev)}$}};")
        lines.append(rf"\draw[->, dashed, mygreen] (i{k}) -- node[rt, text=mygreen, pos=0.32] {{${rate(r)}$}} (c);")
    for k, (r, nxt) in enumerate(outs):
        y = ((len(outs) - 1) / 2 - k) * dy
        lines.append(rf"\node[st] (o{k}) at ({dx},{y:.3f}) {{${state(nxt)}$}};")
        lines.append(rf"\draw[->] (c) -- node[rt, pos=0.68] {{${rate(r)}$}} (o{k});")
    lines.append(r"\end{tikzpicture}")
    return "\n".join(lines) + "\n"


def number(x, digits: int = 4) -> str:
    return f"{x:.{digits}f}"


def macros(values: dict) -> str:
    """\\newcommand definitions for numbers quoted in the slide text."""
    return "".join(f"\\newcommand{{\\{name}}}{{{value}}}\n" for name, value in values.items())
