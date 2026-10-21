# Open Jackson Queueing Network in a Hospital

A reproduction of a group presentation for *Riset Operasi B* (Operations Research B),
Department of Mathematics, Universitas Gadjah Mada. The presentation reviews the paper
*On the Application of the Open Jackson Queuing Network in Hospital* by Arum Helmi Manggala Putri and
Irwan Endrayanto, then builds a simplified three-station model of its own.

## Findings

**The paper's numbers reproduce, after two undocumented steps.**

- Its rate table is the arrival and departure counts divided by 180, i.e. a three-hour observation in
  minutes. The utilisations ρ are truncated, not rounded (0.4722 / 0.7722 = 0.61150 is printed as 0.6114).
- The Moore–Penrose pseudo-inverse printed in the slides is reproduced to within 2.5 × 10⁻⁴.
- The pseudo-inverse solution itself does **not** give the paper's routing probabilities (P₁₂ = 0.4184,
  not 0.4436). They match to within 4 × 10⁻⁷ once the two branches out of stations 1, 2 and 3 are
  rescaled to sum to 1. This step is not stated.

**Two of the paper's headline numbers are wrong.**

- *Mean number of patients, 39.* This is the sum of the per-station means after each was rounded to an
  integer. The unrounded sum is **37.49**.
- *Mean time in the system, 54.61 minutes.* This adds the time per visit of all six stations, as if
  every patient visited every station. Each patient visits four (registration, one examination room,
  one doctor's room, pharmacy). Little's law, W = L / λ₁, gives **36.88 minutes**.

**The paper's routing is inconsistent with its rates.** Feeding the reported routing probabilities into
Jackson's traffic equations gives station 5 a utilisation of **1.031**, so that doctor's room would
never reach a steady state. The routing system is also rank-deficient (rank 8 for 9 unknowns), so
P₂₄, P₂₅, P₃₄ and P₃₅ are not identifiable from the data. The pseudo-inverse only picks the
minimum-norm solution.

**The group's own model needed rederiving.** The slides model registration → doctor (1/3) or
laboratory (2/3), with the laboratory sending 4/5 of its patients on to the doctor. Only one of the
eight balance equations in the slides was correct (the interior case n₁, n₂, n₃ > 0), and six of the
seven hand-drawn transition diagrams had wrong or missing arrows. Here the transitions are generated
by one function, and the equations, diagrams and transition table are all derived from it. The
product-form solution, with λ₂ = 13λ₁/15 and λ₃ = 2λ₁/3, satisfies all eight equations
(checked symbolically with SymPy). The seven original equations that differ fail the same check.

## Results

**Paper (six stations)**

| Quantity | Paper | Reproduction |
|---|---|---|
| λᵢ, μᵢ, ρᵢ (Table 2) | – | identical (counts / 180, ρ truncated) |
| Pseudo-inverse A⁺ | printed to 4 decimals | max. difference 2.5 × 10⁻⁴ |
| Routing probabilities | P₁₂ = 0.4436053, … | A⁺b normalised per station, max. difference 4 × 10⁻⁷ |
| Mean number in system L | 39 | 37.49 (39 after rounding each station) |
| Mean time in system W | 54.61 min | 36.88 min (Little's law); 54.61 is the sum over all six stations |
| ρ₅ under the traffic equations | – | 1.031 (unstable) |

**Three-station model.** The slides give no rates, so the example uses illustrative ones (patients per
minute): λ₁ = 1.0167 and μ₁ = 1.3056 from the paper's registration, μ₂ = 1.2000 from its two doctor's
rooms pooled, and μ₃ = 1.5833 from its two examination rooms pooled.

| Station | λᵢ | ρᵢ | Lᵢ product form | Lᵢ truncated CTMC (≤ 50 per station) | Lᵢ simulation (5 × 200 000 min) |
|---|---|---|---|---|---|
| 1 Registration | 1.0167 | 0.7787 | 3.519 | 3.519 | 3.514 ± 0.035 |
| 2 Doctor | 0.8811 | 0.7343 | 2.763 | 2.763 | 2.780 ± 0.014 |
| 3 Laboratory | 0.6778 | 0.4281 | 0.749 | 0.749 | 0.748 ± 0.003 |
| Total L | | | 7.031 | 7.031 | 7.042 |

W = L / λ₁ = 6.92 min, and the simulation, which times each patient, gives 6.93 ± 0.05 min. Adding the
three stations' time per visit, as the paper does, gives 7.70 min. The truncated chain is solved
without the product form, and its stationary probabilities for states with at most 5 patients per
station agree with it to a relative 3.7 × 10⁻⁶.

All numbers are in [results/paper.json](results/paper.json) and
[results/hospital.json](results/hospital.json).

## Changes to the slides

Corrections:

- **Compilation.** Two stray `\\` line breaks stopped the deck from compiling with TeX Live 2024.
- **Open Jackson definition.** The external rate is λ₀ᵢ ≥ 0 (it said "λ ≥ 0"). The mean service time is
  1/μᵢ (it said E[Wᵢ]). The exit probability pᵢ₀ = 1 − Σⱼ pᵢⱼ is now defined.
- **Jackson's theorem.** The product form had (i − ρᵢ) for (1 − ρᵢ). It now states the stability
  condition ρᵢ < 1.
- **M/M/1.** Arrivals are Poisson and service times exponential (it said arrival and departure *times*
  are Poisson). The local balance is now given for every n (λPₙ₋₁ = μPₙ, so Pₙ = (1 − ρ)ρⁿ) instead
  of only P₁ = (λ/μ)P₀, with L and W added. Fixed the "Daigram" typo.
- **Assumptions.** "No balking" was a copy of the "no reneging" definition. Balking means not joining
  a queue on arrival.
- **State space.** nᵢ counts the patients waiting *and* in service at station i, not only those
  waiting.
- **Transition table.** Added the condition under which each transition can happen.
- **Transition diagrams.** The seven hand-drawn diagrams are replaced by eight generated ones; the case
  (0, 0, 0) was missing. Errors in the originals:
  - (n₁, n₂, n₃): the ⅘μ₃ arrow went to (n₁, n₂, n₃−1) instead of (n₁, n₂+1, n₃−1).
  - (0, n₂, 0): the ⅕μ₃ inflow from (0, n₂, 1) was missing.
  - (0, 0, n₃): the arrow to (1, 0, n₃) was labelled ⅕μ₃ instead of λ₁, and the inflow from
    (0, 0, n₃+1) was labelled μ₃ instead of ⅕μ₃. The μ₂ inflow from (0, 1, n₃) was missing.
  - (n₁, n₂, 0): it had an impossible ⅘μ₃ outflow (the laboratory is empty) instead of
    ⅓μ₁ → (n₁−1, n₂+1, 0). Four inflows were missing.
  - (0, n₂, n₃): both μ₁ inflows came from the wrong states, and the ⅘μ₃ inflow was missing.
  - (n₁, 0, n₃): the ⅔μ₁ outflow was missing, the ⅔μ₁ inflow came from the wrong state, and three
    inflows were missing.
  - (n₁, 0, 0) was correct.
- **Balance equations.** All eight are regenerated. In the originals:
  1. (0, 0, 0): the right side used ⅓μ₁ and ⅔μ₁, which cannot fire when station 1 is empty. It
     should be μ₂P₀₁₀ + ⅕μ₃P₀₀₁.
  2. (n₁, 0, 0): it had no "=", and the λ₁P(n₁−1, 0, 0) term was missing.
  3. (0, n₂, 0): the left side had μ₁ instead of μ₂, the ⅕μ₃ term was missing, and "μ" should be μ₃.
  4. (0, 0, n₃): μ₃ should be ⅕μ₃, and the μ₂ term was missing.
  5. (n₁, n₂, n₃): correct.
  6. (n₁, n₂, 0): the left side had ⅔μ₁ instead of μ₁, and three inflows were missing.
  7. (0, n₂, n₃): both μ₁ terms came from the wrong states, two rates lacked their μ, and the ⅘μ₃
     term was missing.
  8. (n₁, 0, n₃): the left side had ⅓μ₁ instead of μ₁, the μ₁ term came from the wrong state, and
     three inflows were missing.

Additions:

- A third research question, on whether the paper's numbers reproduce.
- A section on the product-form solution, with traffic equations, stability conditions, L and W.
- A numerical example, plus a verification slide that compares the product form with the truncated
  Markov chain.
- A section reproducing the paper's rate table, routing probabilities, L and W.
- The conclusion, which was empty (its commented-out text was the paper's conclusion).
- The table of contents is set in two columns, to fit the nine sections.

The cover no longer lists the members' student ID numbers, here or in `original/`.

## Repository layout

| Path | What it is |
|---|---|
| `src/openjackson/network.py` | Jackson networks: transitions, traffic equations, product form, L and W |
| `src/openjackson/hospital.py` | The slides' three-station model and its symbolic balance equations |
| `src/openjackson/ctmc.py` | Stationary distribution of the truncated Markov chain, without the product form |
| `src/openjackson/simulate.py` | Patient-level discrete-event simulation |
| `src/openjackson/paper.py` | Reproduction of the paper's six-station numbers |
| `src/openjackson/tex.py`, `__main__.py` | Generates `results/` and the slides' `generated/` snippets |
| `slides/` | Corrected deck; `generated/` holds the equations, diagrams and tables written by the code |
| `original/` | The deck as submitted, unchanged |
| `tests/` | Unit tests, including the symbolic check of all eight balance equations |

## Reproduce

```bash
uv venv --python 3.13 && uv pip install -e ".[dev]"
.venv/bin/python -m pytest -q
.venv/bin/python -m openjackson          # about 30 s: writes results/ and slides/generated/
cd slides && latexmk -pdf main.tex
```

## Limitations

- The paper's data come from its tables as shown in the slides. The 180-minute window is inferred
  from them, since every rate equals a count divided by 180.
- The paper estimates each μᵢ from departures counted in the same window. Departures measure
  throughput, not service capacity, and they exceed arrivals at every station. The reproduction keeps
  the paper's method.
- The rates of the three-station example are illustrative, not data.

## Credits

The presentation (Riset Operasi B, Kelompok 4) was made by Ridho Maulana Asrofi, Aprila Adzkiyan Najah,
Nala Kirana Prastiwi, Delia Noor Mavianti and Moh Hasbi Alfarizi. This reproduction is by
Moh Hasbi Alfarizi. The hospital network, its data and the routing method are from Putri and
Endrayanto. The code is MIT-licensed. The slides and their figures belong to their authors.
