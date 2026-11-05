# Initial private release

Prepared November 4, 2026 as a contained snapshot with a fresh repository history. The original private research repository is unchanged.

## Packaging

- All five papers, selected supporting proofs, bound-generation sources and independent checkers.
- A portable small-example runner for counting and both bounds in B and E.
- Large historical certificates, optimization checkpoints, matrix caches, the ledger and unrelated code are omitted.
- Original code is GPLv3 only; original writing and figures are CC BY 4.0, with third-party exclusions recorded separately.

## Manuscript corrections

- A now describes the paired-weight bound on closed walks, qualifies numerical scaling observations, narrows the D-finite statement to the relevant coefficient setting, and removes an unsupported numerical transition location.
- B explicitly separates signed lifting from probability/entropy statements and writes the normalization as `exp(2 H_minus / M_plus)`. Its verification account includes the later audit and distinguishes randomized modular tests from the deterministic construction proof. The upper rational maximum is correctly located in the verification report.
- D replaces an explicitly unverified interpretation of the five-sevenths estimate with a description in terms of the actual retained-rank integrals.
- E spells out the three comparisons needed to transfer either symbol between stacks at fixed parity in the ergodicity argument.
- The shared bibliography points to this selected release and notes that historical appendix references can fall outside its scope.

These edits clarify claims and evidence; they do not constitute a new validation of all mathematical arguments or the record-size certificates. The reproduction guide lists what was executed for this release.
