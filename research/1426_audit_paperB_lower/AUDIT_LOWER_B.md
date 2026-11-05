# 1426 (JB): independent audit of Paper B's lower bound μ ≥ 12.219100857094

**Verdict: CONFIRMED.** On the exported certificate, read through my own code, the final flow J at radius 10 is
positive on all 45,638,378 legal pairs. Every one of the 44,431,506 balance equations (S_10) holds, and so do both
drift equations (D). The entropy satisfies h(J/M) ≥ h_lo = 1.2515001857102438544…, so
e^{2h} ≥ 12.219100857094 = 6109550428547/500000000000, with slack 1.34·10⁻¹⁴ in h. The bound depends on Theorem F
(Paper B §5.3; 1252 Theorem F′), which I did not re-prove (§6 below). Everything else listed in the brief was checked
on 100 percent of the data.

Here μ is the meander growth rate. J(w, ξ) is the flow value of window w and move ξ, β is the reference law, and
x = J/β is the density. ν̄(w) = Σ_ξ x(w, ξ), and M = Σ J is the total mass.

## 1. Independence and inputs

- **What I read.** Paper B (`s2_tape.tex`, `s5_flows.tex`, `s6_lower.tex`), `1421…/SPEC.md`, `1390…/AUDIT.md` and
  its log `logs/lowerB_law.log`, and `1252…/PROOF.md`. I read no project checker or producer code: none of 1019,
  1204, 1228, 1252/1253 scripts or 1421/scripts. I opened no pickle.
- **What I wrote.** All code is my own and lives in this folder: `common.py`, `law.py`, `prep.py`, `parent.py`,
  `walks.py`, `flowcheck.c`, `bound.py`, `primtest.c`, `primtest.py` and `run_all.sh`.
- **Inputs.** `1421_export_lower_certificate/export/*`. Before use, all nine files were verified against that folder's
  `SHA256SUMS`; all are OK (`logs/export_sha256.log`). The two main files are:
  - `parent_r6.bin`: `5afa253a51f8161a81c56aa302dfe1aa60a24fe439588dc081ca30612951832a`;
  - `law.json`: `aaa797ef…c9aae5`.
  The five `dirs_*.tsv` files are as listed in SPEC §11.1.
- **Budget.** At most 4 threads. Peak resident memory was 2.48 GB (`logs/flow_r10_F4.time`). A complete radius-10
  pass, with every pair and every row, takes 8 s wall time. The Python stages take under 1 minute. No file is over
  1 GB, nothing was written outside this folder, and no process is left running.

## 2. Results by task

**1. The law (C1).** `law.py` rebuilds the law from the parameters of Paper B §6.1 in exact rationals.
- The stationary edge law comes from an exact Gaussian solve of the stationarity system. SPEC's closed form is used
  only as a cross-check, and it agrees.
- I used two readings of the perturbation minimum:
  - Paper B's prose: the orbit of the middle heights under reversal and reflection, with off-strip paths counted as 0;
  - SPEC's group labels.
  They give identical tables.
- All 298 entries equal `law.json` exactly (ρ 22, π 42, quartets 80, quintets 154).
- The properties hold: all 154 quintets are positive (minimum 7.520605100194493·10⁻¹¹), the marginals are consistent
  on both sides at every level, the mass is 1, the law is symmetric under reversal and reflection, and the chain on
  4-height states is irreducible.
- This agrees with 1390 (`logs/lowerB_law.log`): the same counts 22/42/80/154, the same minimum to all printed
  digits, and the same properties. 1390 printed no table values, so the comparison is at that level.
- Log: `logs/law.log`.

**2. The parent (C2, C3)**, exact integer arithmetic (`parent.py`, `logs/parent.log`).
- The file has 226,079 records, and D has 14,308 bits.
- The records are exactly the legal pairs at radius 6, in the canonical order (58,450 windows).
- None is zero, and 7,708 are negative.
- Σ N = D, and Σ N Δa = Σ N Δb = 0.
- All 220,110 rows (S_6) have residual exactly 0. Each row is multiplied by D·β(Y), its β ratios are exact rationals,
  and denominators are cleared before an integer comparison.
- The same rows also pass in the C checker modulo the primes (`logs/flow_r6_parent.log`).

**3. The walks (C4)**, in `walks.py` (`logs/walks.log`).
- All 968,894 directions were replayed with the moves of Paper B Table 1: 542,910 + 24,576 + 139,264 + 131,072
  + 131,072.
- For each one, `walks.py` asserted:
  - legality (interior head, move applicable, new word in the strip);
  - closure;
  - protection (first 3 letters, last 3 letters and letter sum unchanged on every tape);
  - every source tape fits the stage radius r_s ≤ 10;
  - zero drift;
  - the coefficient columns agree exactly (hex = num/2^pow ≠ 0).
- 0 failures. The 3,685,260 edges have the walk shapes and counts of `stages.json` in all five stages.
- The C scatter asserts again that every source fits its stage radius, and that every completion is a window at
  which the move is legal.

**4. The final flow (S, D, P)**, in `flowcheck.c`.
- **Construction.**
  - J = J_P↑ + Σ_d c_d Γ_d, as in SPEC §6.2: each edge adds c_d β(A_d)/β(T) to the density at every completion of
    its source tape, and the parent density is lifted unchanged.
  - Densities are exact fixed-point sums. Each term is ⌊term·2⁸⁸⌋, so a pair hit by n terms has
    lo ≤ x·2⁸⁸ ≤ lo + n.
  - At most 35 terms hit any pair. An overflow guard proves that no 128-bit sum can overflow.
- **Changed sets.** The counts per stage equal SPEC §6.3 exactly. Changed pairs / changed windows:

  | stage | pairs | windows |
  |---|---|---|
  | 0 | 3,210,357 | 829,352 |
  | 1 | 621,335 | 372,479 |
  | 2 | 6,165,672 | 3,053,933 |
  | 3 | 10,470,289 | 4,596,747 |
  | 4 | 10,423,865 | 5,020,763 |

- **(P), complete direct pass.** All 45,638,378 values at radius 10 have a rigorous lower bound x_lo > 0. The
  smallest density is 0.0451051240485067645, as in SPEC.
- **(S), complete direct pass.** All 44,431,506 rows (S_10) of the final flow F_4 pass, as do all 11,784,514 rows
  (S_9) of F_1 and all 3,125,954 rows (S_8) of F_0. F_2 and F_3 also pass all 44,431,506 rows at radius 10.
  - The row formula is SPEC 7.1 multiplied by β(Y). I re-derived it from Paper B's (S) (§4 below).
  - Each row is evaluated exactly in Z/p for two random 62-bit primes, p₁ = 2539181208593945057 and
    p₂ = 3795189218594708497.
  - F_4 and F_1 were checked again with a second, fresh pair, 4553340606827994481 and 3936932742539198693
    (`logs/flow_second_primes_F4.log`).
  - 0 failing rows and 0 anomalies, in every run.
- **(D).** Σ βxΔa = Σ βxΔb = 0 mod all four primes, at every stage.
- **Mass.** Σ βx mod p equals the closed form M = 1 + Σ_d c_d k_d β(A_d) (computed exactly in `walks.py`), and the
  closed form lies inside the checker's interval [M_lo, M_hi]. M = 0.96713960459460176936685549729213476873473090440…,
  equal to SPEC's 1253 value to all 58 printed digits. Here k_d is the number of moves of walk d and A_d its anchor.

**5. The bound (C6)**, in `flowcheck.c` and `bound.py`.
- The entropy is evaluated in fixed point, with directed rounding at every step:
  - β enclosed by 64-bit mantissas, relative width at most 4·10⁻¹⁸;
  - my own logarithm: an atanh series in Q64 with a rigorous tail, ln 2 enclosed two ways;
  - the bounds H_lo = Σ β_lo x_lo L_lo(ν̄_lo/x_hi) and M_hi = Σ β_hi ν̄_hi, in Q120.
- Then h_lo = H_lo/M_hi is an exact rational. `bound.py` checks 2h_lo ≥ U, where U is an exact rational upper bound
  of ln(6109550428547/500000000000) (atanh series plus tail). Result: True.
  - h_lo = 1.2515001857102438544530018306544176856208.
  - The project's value is 1.25150018571024385952…; mine is 5.1·10⁻¹⁸ lower, because my rounding is a little looser.
  - The slack is 1.3394·10⁻¹⁴ in h.
  - A Taylor lower bound gives e^{2h_lo} ≥ 12.219100857094327.
- **The primitives were tested against exact values** (`logs/primtest.log`):
  - 20,000 logarithm cases: no case above the true value, largest gap 3.0·10⁻¹⁸;
  - 5,000 random words of 13 to 25 letters: β enclosures valid, and β mod p exact.
- **All five stage values reproduce** (`logs/bound_F*.log`). Each recorded value e^{2h} ≥ … is confirmed, and each
  h_lo lies within 5·10⁻¹⁸ below the recorded h:

  | stage | radius | e^{2h} ≥ | slack in h |
  |---|---|---|---|
  | F_0 | 8 | 12.213777757442 | 1.36·10⁻¹⁴ |
  | F_1 | 9 | 12.214596070132 | 3.73·10⁻¹⁴ |
  | F_2 | 10 | 12.216487956539 | 1.97·10⁻¹⁴ |
  | F_3 | 10 | 12.218343762668 | 1.96·10⁻¹⁴ |
  | F_4 | 10 | 12.219100857094 | 1.34·10⁻¹⁴ |

**6. Damaged data are rejected** (`logs/damage.log`).
- **A flipped walk.**
  - Replay rejects it ("source does not fit" or "not closed").
  - When a non-closed walk (L,R → L,L in stage 0; the last R → L in stage 4) is forced into the flow anyway: 14 rows
    (S_8) fail, and 504 rows (S_10) fail. The drift is nonzero mod both primes.
- **A parent value perturbed by +1/D.**
  - C3 fails (mass and both drifts).
  - The exact C2 check finds 8 failing radius-6 rows.
  - The radius-10 check finds 2,032 failing rows and a nonzero drift.
- **A wrong law entry.**
  - One quintet numerator +1 in a copy of `law.json`: 1 mismatch, and `law.py` exits 1.
  - A law built with ζ = 1/7: all 220,110 exact parent rows fail, all 44,431,506 radius-10 rows fail, and the drift is
    nonzero.
- **Controls.**
  - A coefficient changed so that hex and num/2^pow disagree is rejected.
  - A coefficient multiplied by 2²⁰ (stage 4) or 2⁴ (stage 0) still passes (S), (D) and (P), as the theory says it
    must for any coefficients. The bound then fails: e^{2h} ≥ 2.51, and h short by 2.1·10⁻¹² respectively. So the
    entropy check is live.

## 3. Coverage

| item | coverage | method |
|---|---|---|
| C1 law | 298/298 entries, all properties | exact rationals |
| C2 parent rows | 220,110/220,110 | exact integers |
| C3 parent mass, drift, domain | all 226,079 | exact integers |
| C4 walks | 968,894/968,894 | exact replay |
| (S_10) of F_4, also of F_2 and F_3 | 44,431,506/44,431,506 | exact in Z/p, 2 primes (4 for F_4) |
| (S_9) of F_1 | 11,784,514/11,784,514 | exact in Z/p, 4 primes |
| (S_8) of F_0 | 3,125,954/3,125,954 | exact in Z/p, 2 primes |
| (P) at radius 10 | 45,638,378/45,638,378, direct pass | rigorous fixed point |
| (D), mass | all pairs | Z/p, plus the exact closed form |
| entropy | all 45,638,378 terms | outward rounding, own logarithm |

So 100 percent of every item was covered; nothing was sampled.

## 4. The balance row, re-derived from Paper B

Paper B §5.3 states (S) as Σ_{w,ξ} J(w,ξ) β_{w,ξ}(Y) = ν(Y°) β(Y)/β(Y°). Here Y is a word of 2r+3 letters with
heights, and β_{w,ξ}(Y) is the chance that move ξ from window w, with the far tape drawn from β, shows Y around the
new head. I worked out the four cases directly:
- **L.** The view after the move is two fresh letters followed by all of w, so β_{w,L}(Y) = β(Y)/β(w) and w = Y[2:].
- **R.** By mirror symmetry, w = Y[:n].
- **Split.** The view is w with the head letter replaced by three letters. It is deterministic, so the factor is 1.
- **Merge.** The view comes from V = (two fresh) w (two fresh), with the triple merged, so the factor is β(V)/β(w).
  The two choices of (w_{r−1}, w_{r+1}) give the merge sources.

Dividing by β(Y) gives SPEC 7.1 term by term. I implemented SPEC 7.1 multiplied by β(Y), so the check needs no
modular inverse. The row set is every leveled word of 2r+3 letters in the strip. The counts (220,110; 3,125,954;
11,784,514; 44,431,506) equal SPEC §2.3. In every run there were 0 anomalies, that is, 0 split or merge sources
missing a legal pair.

## 5. What "exact in Z/p" gives

Each row residual r_Y is a fixed rational number, determined by the data. The primes were drawn with
`secrets.randbits` and Miller–Rabin (13 bases, deterministic below 3.3·10²⁴), after the export hashes had been
verified. The program checks that no prime divides D or a numerator or denominator of a law entry (Python's
modular inverse would raise an error otherwise). The 2^pow factors are invertible because the primes are odd. So
reduction mod p is a ring map on the rationals involved, and r_Y ≡ 0 holds exactly when p divides the numerator of
r_Y.

That numerator has fewer than 10⁷ bits. The ingredients are:
- D, 14,308 bits;
- β factors, at most 84 bits for a quartet and 38 bits for a quintet ratio, over at most 25 letters;
- coefficients, 53-bit numerators with exponents at most 80;
- at most about 400 terms per row.

So at most 1.7·10⁵ primes of 62 bits divide it. Each candidate prime is drawn with probability at most
1550/2⁶⁰, since prime gaps below 2⁶⁴ are less than 1550. A nonzero residual therefore passes one prime with
probability at most 2.3·10⁻¹⁰. It passes both primes of a pair with probability at most 5·10⁻²⁰, and all four with
probability at most 3·10⁻³⁹.

Over all 44,431,506 rows, the chance that some nonzero residual escaped is at most 2.3·10⁻¹² with one pair, and at
most 10⁻³⁰ with both pairs (F_4, F_1). This is a randomized exact test, not a deterministic proof of (S).

(S) also has a deterministic route that does not depend on this pass. It uses:
- C2, exact, here;
- C4, here;
- 1252 Propositions A to C (lifting, closed-walk fields, linearity).

I read those propositions, but I did not re-verify them as a referee. The two routes agree.

## 6. Gap between the theorem and what the checker assumed

1. **Theorem F itself** (Paper B §5.3–5.4; full form 1252 Theorem F′) is assumed, together with what it rests on:
   - Theorem 3 (the tape-machine bijection);
   - Lemma 4 (linear completion, in the unrestricted machine);
   - Fekete's lemma.

   The checker verifies its hypotheses on the data:
   - the law is of memory m = 3 with every path of m+2 = 5 heights positive, and the chain is irreducible;
   - r = 10 ≥ m;
   - (S) on every row of the strip;
   - (D);
   - (P), here J > 0 on every legal pair, which is more than J(w, L) > 0.

   It also checks the normalization: M > 0, and Theorem F is applied to J/M with h(J/M) = H/M. I did not re-prove
   Theorem F. It is proved inside the project (805/808, 1252), and Paper B gives only an outline.
2. **Exactness of (S)** in the direct pass is randomized (§5). The deterministic route uses 1252 Propositions A–C.
3. **What came from SPEC and not from Paper B.**
   - The file formats.
   - The explicit move strings. They match Paper B Table 1.
   - The cylinder formula of a memory-three chain. It is standard, and its marginal identities were verified.
   - The explicit field formula. Paper B §6.2 states it in prose ("add z β(T₀)/β(T) to the density … at every window
     whose central letters … spell T"), and it matches.
   - Evaluating each direction at its stage radius and lifting the density. This is a combinatorial identity: when
     the tape fits r_s, a radius-10 window completes it exactly when its radius-r_s center does.
   - SPEC's law formulas coincide with my Paper B reading (§2, item 1).
4. **Interpretive points that did not matter.**
   - Head cycles are written L,R rather than R,L (the same cycle).
   - The `memory_eta` naming slip (I used η = 3/16 and ζ = 1/8 as Paper B states; the tables match `law.json`).
   - Off-strip paths in the minimum (both readings give the same tables).
5. **Not done.** I did not referee Theorem F or Propositions A–D, and nothing is formalized. The certificate's
   coefficients were taken as given; they play no role in validity. I did not compare against the producer's saved
   changed sets or float increments, since those are inside the npz files, which I did not open. My changed-set
   counts equal SPEC's numbers instead.

## 7. Files and digests

- **Scripts.**
  - Python: `common.py`, `law.py`, `prep.py`, `parent.py`, `walks.py`, `bound.py`, `primtest.py`.
  - C: `flowcheck.c`, `primtest.c`.
  - `run_all.sh`.
- **Logs** (`logs/`):
  - `export_sha256.log`, `law.log`, `parent.log`, `walks.log`, `primtest.log`;
  - `flow_r6_parent.log`, `flow_r8_F0.log`, `flow_r9_F1.log`, `flow_r10_F2.log`, `flow_r10_F3.log`,
    `flow_r10_F4.log` (plus `.time` files);
  - `bound_F0.log` to `bound_F4.log`;
  - `flow_second_primes_F4.log`, `flow_second_primes_F1.log`, `bound_second_primes_F4.log`;
  - `damage.log`, `damage_parent_r6.log`, `flow_r10_bad_scale*.log`.
- **Intermediate binaries.** They were deleted after use and can be rebuilt in about 1 minute.
  `logs/intermediate_sha256.txt` records their SHA-256: `edges_c.bin` `d5726597…af1fe` (identical in two
  independent runs) and `parent_c.bin` `afeccf74…8fe1`.
- **Digests.** `SHA256SUMS` in this folder covers all scripts, logs and small work files.
