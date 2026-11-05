# 1426 FINAL: independent audit of Paper B's lower bound μ ≥ 12.219100857094 (prefix JB)

Auditor folder 1426, written 2026-10-26, heavy track (at most 4 threads; peak 2.48 GB; a complete radius-10 pass
takes 8 s). The full report is in `AUDIT_LOWER_B.md`.

## Attempted

The goal was to confirm or refute Theorem F's application in Paper B §5–6 with my own checker. I wrote it from
Paper B and `1421…/SPEC.md` alone. I read no project checker or producer code, and I opened no pickle. The checker
ran on the hash-verified export `1421_export_lower_certificate/export/*`.

## Verdict: CONFIRMED (certified (exact computation), resting on Theorem F)

- **Law (C1).** I rebuilt it from Paper B §6.1 in exact rationals, under both readings of "smallest probability".
  Both readings give the same tables, and all 298 `law.json` entries agree exactly. All 154 five-paths are positive
  (minimum 7.520605100194493·10⁻¹¹), the marginals are consistent, both symmetries hold, and the chain is
  irreducible. This is the same as 1390.
- **Parent (C2, C3).** Exact integer arithmetic:
  - all 220,110 rows (S_6) have residual 0;
  - Σ N = D, and both drifts are 0;
  - the 226,079 records are exactly the legal pairs at radius 6, in canonical order;
  - 7,708 values are negative.
- **Walks (C4).** All 968,894 directions are legal, closed, protected and drift free, and each fits its stage radius.
  0 failures, and the shapes equal `stages.json`.
- **Final flow.**
  - **(P).** A complete direct pass: all 45,638,378 values at radius 10 are positive. The smallest density is 0.0451051.
  - **(S).** A complete direct pass of every equation:
    - all 44,431,506 rows (S_10) of the final flow F_4 (and of F_2 and F_3);
    - all 11,784,514 rows (S_9) of F_1;
    - all 3,125,954 rows (S_8) of F_0.

    The rows were checked exactly in Z/p for randomly drawn 62-bit primes: two primes, and four for F_4 and F_1.
    0 rows failed. The chance that a nonzero residual escaped is at most 10⁻³⁰ for F_4.
  - **(D).** Both drift sums are 0 mod every prime.
  - **Mass.** The mass equals the exact closed form 1 + Σ c_d k_d β(A_d) = 0.967139604594601769366855497292…
  - **Changed sets.** The counts per stage equal SPEC exactly.
- **Bound.**
  - h_lo = 1.2515001857102438544…, from my own outward-rounded logarithm and fixed-point arithmetic.
  - 2h_lo ≥ an exact rational upper bound of ln(6109550428547/500000000000). Slack 1.34·10⁻¹⁴ in h.
  - All five stage values (12.213777757442 … 12.219100857094) reproduce, each within 5·10⁻¹⁸ of the project's h.
- **Damage tests.** Each damaged input is rejected:
  - flipped or non-closed walks;
  - a parent numerator +1;
  - a wrong `law.json` entry;
  - a wrong ζ;
  - inconsistent coefficient columns.

  Controls behave as the theory predicts: rescaled coefficients keep (S), (D) and (P), but the bound then fails.
- **Closed here.** The SPEC §10 gaps "(S_9)/(S_10) only sampled (0.9 %, 0.73 %)" and "positivity never checked by a
  direct pass over all 45,638,378 values" are now closed. Both checks are complete.

## Depends on

- Theorem F (Paper B §5.3; 1252 Theorem F′). It rests on Theorem 3, Lemma 4 and Fekete's lemma. I assumed these and
  did not referee them.
- Exactness of (S) has two routes:
  - deterministically, it follows from C2 and C4 (both checked exactly here) and 1252 Propositions A–C (proved in the
    project; not refereed by me);
  - the direct pass is a randomized exact test.

## Failed or not done

- No referee pass on Theorem F or on Propositions A–D.
- No deterministic exhaustive exact rational check of (S_10). It would take hours in pure big-integer Python; the
  modular pass replaces it.
- The producer's saved changed sets were not compared, because they are inside npz files, which are not allowed. My
  enumerated counts equal SPEC instead.

## Sharpest open statement and next step

The certificate is now checked in full, apart from Theorem F itself. The next step is an independent referee of
Theorem F′ and Propositions A–D (1252). In particular, the counting step in the unrestricted machine needs checking,
with the strip used only for the stationary middle part.

## What this says about q > 0 and multi-loop systems

Nothing. This audit concerns the single-loop bound only.

## Ledger rows

| ID | Finding | Status | Best folder | Other folders | Goals | Class | Route | Object |
|---|---|---|---|---|---|---|---|---|
| JB1 | Independent confirmation of Paper B's lower bound mu >= 12.219100857094 = 6109550428547/500000000000: a new checker written from Paper B and the 1421 SPEC only (no project code, no pickle), run on the hash-verified export, rebuilds the law (298 entries equal law.json exactly, 154 positive five-paths, minimum 7.520605100194493e-11), checks the radius-6 parent exactly (all 220,110 rows residual 0, sum N = D, drifts 0, 7,708 negative values), replays all 968,894 walks (closed, legal, protected, fitting, drift free; 0 failures), and bounds the entropy of the final radius-10 flow with its own outward-rounded logarithm: h >= 1.2515001857102438544 (slack 1.34e-14 in h; all five stage values 12.213777757442 to 12.219100857094 reproduced within 5e-18 in h); rests on Theorem F (Paper B 5.3, 1252 Theorem F'), not refereed here | certified (exact computation) | 1426_audit_paperB_lower | 1421_export_lower_certificate, 1228_combined_lower, 1253_referee_lower_computation, 1252_referee_lower_theory, 1390_cross_model_audit | mu_low | standing | grammar | all |
| JB2 | First complete direct checks of the final lower-bound flow, closing the gaps of 1421 SPEC section 10: positivity on all 45,638,378 radius-10 pairs by a direct pass (smallest density 0.0451051240485067645), and every balance equation one by one (all 44,431,506 rows of (S_10) for the final flow and for stages 2 and 3, all 11,784,514 rows of (S_9) for stage 1, all 3,125,954 rows of (S_8) for stage 0) exactly in Z/p for random 62-bit primes (four primes for the final flow; escape probability at most 1e-30), plus both drifts 0 mod p and the enumerated mass equal to the exact closed form 1 + sum c_d k_d beta(A_d) = 0.96713960459460176936...; the changed sets per stage equal SPEC (621,335; 6,165,672; 10,470,289; 10,423,865 pairs); deterministic exactness of (S) still rests on the parent check plus 1252 Propositions A-C | certified (exact computation) | 1426_audit_paperB_lower | 1421_export_lower_certificate, 1253_referee_lower_computation | mu_low | standing | grammar | all |
| JB3 | Independent checker for the Paper B lower-bound certificate (Python exact parent and walk replay; C pass over all pairs and rows with fixed-point densities, directed-rounding mantissas, an own atanh logarithm and Montgomery arithmetic mod two random primes; 8 s and 2.5 GB at radius 10), with damage tests: a non-closed walk forced into the flow makes 14 (S_8) or 504 (S_10) rows fail, a parent numerator +1 makes 8 exact radius-6 rows and 2,032 radius-10 rows fail, a wrong law.json entry or zeta = 1/7 is rejected (all rows fail), and rescaled coefficients keep (S), (D), (P) while the bound fails | infrastructure | 1426_audit_paperB_lower | 1421_export_lower_certificate | mu_low | infra | grammar | all |
