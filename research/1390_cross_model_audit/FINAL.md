# 1390: cross-model audit of the headline certified bounds (prefix HP)

**Task.** HANDOFF_CODEX_2026-10-26.md, Task 1. Confirm or refute four headline bounds with checkers written
independently, from the papers' theorem statements and format descriptions only. I did not read or run any project
checker or producer code. Heavy track: at most 4 threads; peak 2.6 GB; the longest run took 25 s. Details are in
`AUDIT.md`; digests in `SHA256SUMS`; logs in `logs/`.

## Results

| bound | statement used | verdict |
|---|---|---|
| Paper E: μ² ≤ 10.157042906664444025 | Lemma 5, Lemma 6 (truncation); format from 602 PROOF.md §6 | CONFIRMED |
| Paper E: μ² ≥ 10.086918840819 | Theorem L, (S) checked directly; Proposition 7 used only to build J | CONFIRMED |
| Paper B: μ ≤ 12.2970 (R_+² = 12.29693385352921289818…) | Lemma 6, Corollary 7, Lemma 8, Proposition 9 | CONFIRMED |
| Paper B: μ ≥ 12.2191 | Theorem F with the §6 constructions | NOT CHECKED (undocumented formats, pickled parent) |

**Status words.** The three confirmations are `certified (exact computation)` by a checker independent of the
producers. Each rests on the theorem named; I did not re-referee the proofs of those theorems. Two points need
care:
- For bound 3, the JSON encoding is undocumented and was inferred. The inference is exact: the stored R_+ is
  reproduced as an exact fraction at the stored argmax. By Proposition 9 any positive reading would give a valid
  bound.
- For bound 2, the logarithms are enclosed by my own integer series, with all rounding in the safe direction. No
  libm or mpmath is involved.

## Main numbers

- **E upper.** Every one of the 536,805,378 states with stacks of at most 13 symbols passes, exactly. Lowering P by 1
  gives 4 violations. The largest quotient is 3.187011594348 at U = 0, L = 1111000000, c = 0.
- **E lower.** All 8,388,608 balance equations hold exactly, and both drifts are 0. h ≥ 1.155619710060893, so
  e^{2h} ≥ 10.0869188408196; the target …819 is met and …820 is not.
- **B upper.** Of the 2,097,152 formal rows, the maximum equals the stored R_+ exactly, at row 697690. In 2,000
  actual tapes of up to 150 letters, no actual row value exceeds its formal row value.
- **Damaged certificates.** Every damaged certificate or wrong input was rejected: 27 of 27 (9 for the E upper
  bound, 10 for the E lower, 8 for the B upper).

## Paper-level issues found (none invalidates a bound)

1. Paper E gives no byte format for the upper certificate (only 602 PROOF.md §6 does).
2. Paper E records no hash for `flow_r10.bin`; mine is `20bcd572…19fd39e6`.
3. Paper E, Theorem L, Step 2: the shift invariance at a fixed parity needs a three-move composition that is not
   written. I checked that it holds.
4. Paper B, upper: the JSON encoding is undocumented, and R_+ is not in the certificate file, contrary to §7.
   Whether the padding letters belong to the gram sides is not stated (harmless by Lemma 8).
5. Paper B, lower: the reference-law text needs two interpretive choices, listed in AUDIT.md §4. Under them, the law
   rebuilt from the text satisfies all hypotheses of Theorem F: 154 positive five-paths, consistent marginals,
   irreducible. The lifting, closed-walk and sparse-refinement constructions are prose, not numbered results. The
   radius 9 and 10 balance equations are only sampled anywhere in the project.
6. Metadata naming slip: `"memory_eta": "1/8"` in `combined_cycle_r8.npz`, while η = 3/16.

## Open / next step

Bound 4 needs three things before an outside-style check is possible:
- a documented, non-pickle export of the radius-6 parent and of the walk descriptors;
- numbered statements of 1252 Propositions B–D in Paper B;
- preferably, the final radius-10 flow values, so that (S) can be checked exhaustively in exact arithmetic.

## What this says about q > 0 and multi-loop systems

Nothing.

## Ledger rows

| ID | Finding | Status | Best folder | Other folders | Goals | Class | Route | Object |
|---|---|---|---|---|---|---|---|---|
| HP1 | Cross-model independent exact check (own C code from Paper E Lemmas 5-6 and the BICCERT1 format only) of the radius-12 bicubic upper certificate cert_r12.bin (sha256 3e71e227...a219, matches 602): all 536,805,378 short states pass, 4 fail with P-1, largest quotient 3.187011594348 at U=0, L=1111000000, c=0, so z_n <= (3187011595/10^9)^(2n) and mu^2 <= 10.157042906664444025; truncation lemma confirmed on 2e7 random long states; 9 damaged copies rejected | certified (exact computation) | 1390_cross_model_audit | 602_bicubic_upper_bound, 1367_bicubic_note | relatives | standing | window_weight | model |
| HP2 | Cross-model independent exact check (own code from Paper E Theorem L, Proposition 7 and the BICFLOW1 format only) of flow_r10.bin (sha256 20bcd572...19fd39e6): (P), (D) and all 8,388,608 balance equations (S) hold exactly when checked directly by applying each move to concrete stacks with two hidden reference symbols per stack; h >= 1.155619710060893 by a downward-rounded integer atanh series, so mu^2 >= 10.0869188408196 >= 10.086918840819 (target 10.086918840820 fails); radii 4 and 8 reproduce 10.030980878326 and 10.079276064125; damaged flows rejected | certified (exact computation) | 1390_cross_model_audit | 1367_bicubic_note, 1377_referee_bicubic_lower | relatives | standing | entropy_flow | model |
| HP3 | Cross-model independent check (own code from Paper B Lemma 6, Corollary 7, Lemma 8, Proposition 9 only) of the radius-9 paired certificate r9_run2_cand.json (sha256 332b0346...2246): float filter with rigorous error bound over all 2,097,152 formal rows plus exact Fraction evaluation of the 4 rows within 1e-12 gives maximum exactly the stored R_+ at row 697690, R_+^2 = 12.29693385352921289818... <= 12.2970; two separate implementations agree to 2.6e-16 on 3,000 rows; on 2,000 actual tapes up to 150 letters every actual row value is at most its formal row; damaged tables rejected; the JSON encoding is undocumented and was inferred (any positive reading is valid by Proposition 9) | certified (exact computation) | 1390_cross_model_audit | 1201_upper_structures, 1251_referee_upper_12297 | mu_up | standing | window_weight | all |
| HP4 | The Paper B lower bound 12.2191 was not checked independently: its certificate is a pickled radius-6 parent plus 968,894 int64-encoded walk descriptors with no format description outside project code; what was checked: the memory-three reference law rebuilt from Paper B Section 6.1 alone (reading: off-strip paths count 0 in the minimum; orbit of the middle heights under reversal and reflection) has 154 positive five-paths, consistent marginals, both symmetries and an irreducible chain, so the law hypotheses of Theorem F hold for the law as described; open: a documented export of the parent and descriptors and numbered statements of the lifting and closed-walk constructions | infrastructure | 1390_cross_model_audit | 1228_combined_lower, 1253_referee_lower_computation, 1252_referee_lower_theory | mu_low | infra | grammar | all |
