# Meandric research

Papers, selected proofs, and programs for studying the growth of meander numbers. This is a contained research snapshot by Carl Johan Ragnarsson, prepared on November 4, 2026. The repository is initially private.

**[Read the papers (PDFs)](papers/) · [Start with Paper B: bounds](papers/B_tape_bounds.pdf)**

The main contribution is the tape method in Paper B: encode a meander as a closed sequence of local edits, obtain upper bounds from paired weight inequalities, and obtain lower bounds from the entropy of stationary flows. The paper reports **12.219 ≤ μ ≤ 12.297**, where μ is the exponential growth constant of closed meanders counted by order. Bounds here are rounded outward to three decimal places; Paper B gives the full precision.

| Source | Lower bound | Upper bound | Publication |
|---|---:|---:|---|
| [Albert and Paterson (2005)](https://doi.org/10.1016/j.jcta.2005.02.006) | 11.380 | 12.901 | Journal article |
| [Fox (2026)](https://github.com/mnfoxx/meander-bounds) | 11.608 | 12.632 | Manuscript and code |
| [Caroppo, Da Lozzo, Di Battista and Vihrovs (2026)](https://arxiv.org/abs/2610.11944) | — | 12.806 | Preprint |
| **[This work: Paper B](papers/B_tape_bounds.pdf)** | **12.219** | **12.297** | **Research draft** |

This comparison follows Table 1 of Paper B. The numerical estimate μ ≈ 12.263 is not a proved bound.

This repository contains code to generate and independently check smaller examples of those methods. It does not distribute the large historical certificates or claim a tested reconstruction of every record endpoint. See [Reproducing the bounds](docs/REPRODUCING.md) for the distinction and the current limits.

## Papers

| Paper | Topic | Read | Source |
|---|---|---|---|
| A | Intuition and the open meander problem | [PDF](papers/A_intuition.pdf) | [LaTeX](research/paper_series/A_intuition/main.tex) |
| B | Tape encoding and growth-constant bounds | [PDF](papers/B_tape_bounds.pdf) | [LaTeX](research/paper_series/B_tape_bounds/main.tex) |
| C | Subexponential estimates and bounded-label families | [PDF](papers/C_subexponential.pdf) | [LaTeX](research/paper_series/C_subexponential/main.tex) |
| D | Continuum models and conditional connections to finite meanders | [PDF](papers/D_continuum.pdf) | [LaTeX](research/paper_series/D_continuum/main.tex) |
| E | Bounds for Hamiltonian cycles on bicubic maps | [PDF](papers/E_bicubic.pdf) | [LaTeX](research/paper_series/E_bicubic/main.tex) |

These are research drafts. The full meander asymptotic formula and its predicted critical exponent remain unproved. No external peer review is claimed. C and D have selected supporting proofs here, not their entire historical dependency chains.

## Run the small examples

Use Python 3.12 and a C/C++17 compiler for the tested setup:

```sh
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python scripts/reproduce.py all
```

The examples count small objects, generate fresh upper/lower certificates for B and E, and run separate checkers. Outputs go to `outputs/`, which is ignored by Git. Each example can also be selected individually: `counts`, `b-upper`, `b-lower`, `e-upper`, or `e-lower`.

The small runs give weaker bounds than the headline results. Their purpose is to make the mathematical constructions executable from a clean checkout. Numerical optimization alone is never presented as a certified bound.

## Read the methods

- [Reproduction guide and measured small-run results](docs/REPRODUCING.md)
- [Map of included proofs and code](docs/CONTENTS.md)
- [Changes made for this release](docs/RELEASE_NOTES.md)
- [Attribution and third-party material](docs/ATTRIBUTION.md)
- [Final file checksums](MANIFEST.tsv), checked by `python scripts/check_manifest.py`

Source paths retain their original research folder names so that the paper appendices remain useful. A reference to another historical folder does not imply that folder is included. Large certificates, optimization checkpoints, generated matrices, the research ledger, and unrelated experiments are excluded.

## Licenses

Original code is **GPL-3.0-only**. Original papers, proofs, figures, and documentation are **CC BY 4.0**. Third-party material retains its own terms. See [LICENSE](LICENSE), [GPLv3](LICENSES/GPL-3.0-only.txt), [CC BY 4.0](LICENSES/CC-BY-4.0.txt), and the attribution note.
