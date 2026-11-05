# Reproducing the bounds

The executable examples construct certificates locally and check the resulting inequalities independently. No large input certificates are downloaded or included. The default runs use small radii and do not reproduce the strongest endpoints reported in the papers.

## Tested setup and results

The release examples were run on macOS with Python 3.12.15, NumPy 2.5.4, SciPy 1.18.1, mpmath 1.4.1, and the system C/C++ compiler. Dependency versions are pinned in `requirements.txt`. Linux shared-library naming is supported, but a Linux run has not been performed for this release.

| Command after `python scripts/reproduce.py` | Generation and independent check | Result in the release check |
|---|---|---|
| `counts` | Noncrossing matchings versus the tape machine, orders 1–6 | 1, 2, 8, 42, 262, 1828; also checks the two bounded-label counts and bicubic model |
| `b-upper` | Paired weights at radius 3, starting from zero log-weights; exact rational evaluation of all 512 formal rows | μ ≤ approximately 12.442561071797494; exact rational upper bound written to `outputs/b_upper_bound.json` |
| `b-lower` | Radius-3 stationary-flow optimization, exact rational reconstruction, then a separate balance/drift/entropy checker | μ ≥ 11.843448806 |
| `e-upper` | Radius-3 two-stack weight solve, integer certificate, then a separately implemented full-state checker | μᵦ² ≤ 10.592570872129 |
| `e-lower` | Radius-4 reversible-flow iteration, integer certificate, then exact balance and directed-rounding entropy checks | μᵦ² ≥ 10.0309808783260 |

Optimizer output and last digits may vary with libraries and machine arithmetic. The exact or outward-rounded check of each newly generated candidate is authoritative. B's μ is per order; E's μᵦ is per vertex, so the table reports its square.

The generated pickle in the B lower example is loaded only after it is produced locally by the preceding step. Historical checkpoint formats and executable inputs are not downloaded by the example runner.

## Paper B upper bound

The included producer is `research/1201_upper_structures/run_big.py`, with the sparse model, lifting, quantization and native evaluation dependencies. The small wrapper compiles the native helper, creates an initial vector, searches for weights, and writes a candidate JSON. `research/1251_referee_upper_12297/run_formal.py` independently evaluates every row in rational arithmetic.

The resulting exact maximum R gives μ ≤ R² by the paired-weight argument. The stored `r9_run2_verify.json` is a historical report, not its 78 MB weight certificate.

The paper's radius-9 endpoint came from continuation of earlier optimized weights and a time-limited run. This subset supplies the method and the model, but does not contain the historical warm starts or a validated deterministic recipe that reproduces those final digits. Increasing the radius requires a larger search and a fresh exact check.

## Paper B lower bound

The small complete route uses `research/805_entropy_certificates/exact_cert.py` and its independent verifier. It solves a stationary-flow problem, reconstructs rational values, checks balance and zero drift exactly, and bounds entropy with interval arithmetic. Its weaker bound illustrates the same entropy principle.

Selected advanced sources include the memory-three reference law, signed-parent reconstruction, protected-cycle directions, sparse refinements, template search and portable export. The main theorem and construction identities are in `research/1252_referee_lower_theory/PROOF.md`; the export specification is in `research/1421_export_lower_certificate/SPEC.md`.

The exact record 12.219100857094 uses an optimized radius-6 parent, a completed radius-8 flow, and four larger-radius increments. Reconstructing this historical chain needs earlier base-parent and basis data plus the optimization sequence. Those data are excluded here. The advanced source selection is therefore a reference for extending the method, not an independently runnable record-reconstruction command.

The later audit in `research/1426_audit_paperB_lower` can check a full locally generated export. Its balance proof uses exact parent equations and checked cycle identities; its exhaustive modular checks are additional randomized consistency checks. It should not be described as deterministic rational evaluation of every final balance equation.

## Paper E

`research/602_bicubic_upper_bound/bic.cpp` contains both solve and certificate modes. The wrapper first writes a floating-point table, converts it to an integer certificate, and passes that certificate to the independent checker in `research/1390_cross_model_audit/upperE`.

The lower generator consists of `make_cert.py`, `fixpoint.py`, and `flow_opt.py` in `research/1367_bicubic_note`. It starts from a simple positive flow. The wrapper runs 400 iterations at radius 4 and invokes the independent checker in `research/1390_cross_model_audit/lowerE`, including all 2,048 extended-window balance equations.

Larger radii use the same constructions but consume more memory and time. The radius-12 upper and radius-10 lower certificates behind the reported record bounds are not included. A fresh large run must be checked before its result can be called certified.

## Build the papers

Install Tectonic separately, then run `tectonic main.tex` inside each paper directory. The figures and bibliography are included. The checked-in PDFs are built from the release sources. The paper appendices refer to both included and omitted historical folders; consult `CONTENTS.md` rather than treating the release as the complete private research archive.

## Scope of validation

All five small examples were run for this release. Source syntax, native builds used by those examples, PDF compilation and layout, and file checksums were checked. No large certificate verification, record-size optimization, or full proof re-review was performed.
