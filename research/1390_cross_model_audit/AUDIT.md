# 1390: independent audit of the four headline certified bounds (prefix HP)

Auditor: a model of a different family from the producers and earlier referees. I wrote every checker in this
folder from the theorem statements in the papers. For the formats I used only format descriptions: 602 PROOF.md
§5–6, Paper E's appendix, and the JSON keys themselves. I did not read or run any project checker or producer code.
I opened no pickle. All runs used at most 4 threads; peak memory was 2.6 GB.

| # | bound | verdict | evidence |
|---|---|---|---|
| 1 | Paper E upper, μ² ≤ 10.157042906664444025 | **CONFIRMED** | all 536,805,378 short states checked exactly; 0 violations |
| 2 | Paper E lower, μ² ≥ 10.086918840819 | **CONFIRMED** | (P), (D) and all 8,388,608 equations (S) exact; h bounded with rounding in the safe direction |
| 3 | Paper B upper, μ ≤ 12.2970 | **CONFIRMED** | all 2,097,152 formal rows; the exact maximum equals the stored R_+ |
| 4 | Paper B lower, μ ≥ 12.2191 | **NOT CHECKED** | formats undocumented (see §4); only the reference law was rebuilt and tested |

Digests: `SHA256SUMS` in this folder covers the certificates and my checkers. Logs are in `logs/`.

## 1. Paper E, upper bound (Lemma 5 and Lemma 6, `s3_upper.tex`)

- **Input.** `602_bicubic_upper_bound/cert/cert_r12.bin` (a Lexar symlink), 268,435,540 bytes. SHA-256
  `3e71e227…bfca219`, equal to 602's `cert/SHA256SUMS`. Header: r = 12, a = 1783/2000, p = 137/200, t = 97/125,
  R = 3187011595/10⁹, count 4¹³.
- **Checker** (`upperE/hp_upperE.c`). It enumerates every state (U, L, c) whose stacks have at most r+1 = 13
  symbols. For each, it applies the four moves literally and recounts E and T on the whole stacks. It reads the views
  off the stacks and the parity-1 entries through the 0↔1 swap. Then it tests Σ_y v(y) ≤ R v(x) in unsigned 128-bit
  integers after clearing denominators; a start-up bound shows every value stays below 2¹²⁶.
- **Result** (`logs/upperE_r12.log`, 8 s wall).
  - 536,805,378 states, 0 violations.
  - With P − 1: 4 violations.
  - Largest quotient 3.187011594348 at U = 0, L = 1111000000, c = 0. The tightest state with two full views is
    3.187011581921. Both agree with the paper and with 1367.
  - Every one of the 67,092,481 used entries is positive.
  - R² = 406281716266577761/4·10¹⁶ = 10.157042906664444025 exactly.
- **Second, separate check** (`upperE/hp_upperE_literal.py`, Python Fractions).
  - It computes v from scratch, not as a ratio, on 60,000 random states: long stacks up to 30, states reached by
    random walks from e, and near-tight families. 0 violations; largest quotient seen 3.187011594348.
  - A transfer count gives z₁…z₁₀ = 2, 8, …, 29114128, as in the paper, and z_n ≤ R^{2n} holds for them.
  - Lemma 6 (truncation) was tested on 2·10⁷ random states with stacks up to 39 symbols: the scaled sums match
    their truncations exactly (`logs/upperE_r12_tests.log`).
- **Damage tests** (`logs/upperE_r12_tests.log`, `logs/upperE_r12_filetests.log`). Every damaged copy is rejected:
  - R lowered by 10⁻⁹ in the header;
  - the tightest entry −0.1 %;
  - the tightest two-full-view entry −0.1 %;
  - a successor's entry +0.1 %;
  - a zero entry;
  - truncation by 4 bytes, a trailing byte, wrong magic, a wrong radius field.
- **Radius 7.** It reproduces 3.191740267742 at U empty, L = 000000.

**Gaps (paper level, none affects validity).**
- Paper E does not give the byte format of the upper certificate (it gives the lower one). The format lives only in
  602 PROOF.md §6, which Paper E does not cite for this purpose.
- Lemma 5 requires the inequality "for every state"; the certificate meets that, which is more than the reachable
  states (charge rule) need.
- Nothing had to be assumed beyond the paper and 602 §6.

## 2. Paper E, lower bound (Theorem L and Proposition 7, `s4_lower.tex`; format in the appendix)

- **Input.** `1367_bicubic_note/cert/flow_r10.bin`, 16,777,260 bytes, SHA-256 `20bcd572…19fd39e6`. Neither Paper E
  nor 1367's FINAL records a hash for this file; mine is the first recorded. Header: r = 10, s = 41/50,
  n_A = n_B = 4¹⁰. Every X is positive, and max X = 2⁶⁰ as stated.
- **Checker** (`lowerE/hp_lowerE.py`).
  - It builds J from Proposition 7 in exact integers (the common factor 8·50¹⁹ is dropped).
  - It checks (P) and (D) exactly.
  - It checks (S) without using Proposition 7. For each window, allowed move and each value of the 4 hidden symbols
    (2 below each stack, drawn from the reference chain), it applies the move to the concrete top r+2 symbols. It then
    compares the mass of every extended window after the move with its mass before, exactly.
  - It bounds h from below with an integer fixed-point atanh series, every rounding downward (2⁻⁹⁶ fixed point;
    ln 2 enclosed by partial sums plus a geometric tail). It compares 2h with an upper bound of ln(target) in exact
    rationals.
- **Result** (`logs/lowerE_r10.log`, 25 s, 2.6 GB).
  - 2,097,152 windows and 6,291,456 flow values; (P) holds; the drift sums are exactly 0.
  - All 8,388,608 equations (S) hold exactly.
  - h ≥ 1.155619710060893, so e^{2h} ≥ 10.0869188408196. The target 10.086918840819 is met with slack 3.2·10⁻¹⁴
    in h; the target …820 is not met.
  - Radius 4 reproduces 10.030980878326 (slack 4.6·10⁻¹⁵ in h) and radius 8 reproduces 10.079276064125.
- **Damage tests** (`logs/lowerE_r10_tests.log`).
  - Rejected because (S) fails: pushes without q (all 8,388,608 equations fail and (D) fails); one flow value +1
    (6 equations fail, (D) off by 1); a pop that reveals the wrong fresh symbols (4,194,304 fail).
  - Rejected at load: a zero X, a truncated file, a trailing byte, wrong magic, a wrong radius field.
  - Exchanging the A and B blocks still satisfies (S) and (D), as Proposition 7 predicts, but gives only h ≥ 0.856,
    so the target is not met. This confirms that my reading of the index formula is the producer's.

**Gaps (paper level).**
- The format does not say whether the 64-bit X are signed. This does not matter, since every X < 2⁶¹.
- Step 2 of Theorem L derives the shift invariance from "a push on U with a push on L at the two parities". Moving a
  symbol equal to the current parity needs a three-move composition: push, cross at the other parity, then compare
  with a push on L. I checked that this composition gives f(0U,L,0) = f(U,0L,0). The argument is correct but terse.
- No hypothesis is missing. Rounding direction is not an issue: everything is exact except the logarithms, which I
  rounded in the safe direction myself.

## 3. Paper B, upper bound (Lemma 6, Corollary 7, Lemma 8, Proposition 9; Table 4)

- **Input.** `1201_upper_structures/runs/r9_run2_cand.json` (a Lexar symlink), SHA-256 `332b0346…dcee52246`, equal
  to Paper B's appendix and 1201's SHA256SUMS. The claimed R is in `r9_run2_verify.json`.
- **Format.** No document describes the JSON encoding. I inferred it from the keys and the symmetries Paper B states:
  - letter + = 1;
  - window, cut and edit-window codes are most-significant-bit first, leftmost letter first;
  - gram codes are most-significant-bit first in reading order (from the far end toward the head);
  - `field_lookup[cut] = [index, exponent]`, with E = E[index]^exponent and exponent 0 meaning 1;
  - `edit_lookup` key = window·8 + triple.

  The stored tables then show exactly the stated symmetries (0 mismatches). Two facts make this safe:
  - By Proposition 9, any reading that gives positive tables with reciprocal factors gives a valid bound. A wrong
    reading could only change the number.
  - The other readings give 10.25 or 6.07 (`logs/upperB_readings.log`); the bit order of G is irrelevant by symmetry.
- **Checkers.** Two implementations, written separately.
  - `upperB/hp_rowsB.c`, a float filter. It makes every move and every predecessor literally on a padded word and
    computes W and Z from the whole words, cancelling equal gram factors exactly. Rigorous error bound below
    10⁻¹⁴ against a margin of 10⁻¹².
  - `upperB/hp_upperB.py`, an exact Fraction evaluator.
- **Result** (`logs/upperB_exact.log`).
  - All 2,097,152 rows were filtered; 4 rows lie within 10⁻¹² of R. Exact evaluation of those 4: the largest equals
    the stored R exactly, at row 697690 = `-+-+-+-+--[+]-+-+-++-+-`, as stored. No row exceeds it.
  - R_+ = 3.5066984263733334…; R_+² = 12.29693385352921289818… (`logs/upperB_R2.log`) ≤ 12.2970.
  - The two implementations agree to 2.6·10⁻¹⁶ relative on 3,000 random rows.
- **Direct test of the reduction** (`logs/upperB_tapes_big.log`). On 2,000 actual tapes of up to 150 letters
  (random walks from the root), I used full-tape weights and only the moves the machine really has. Every actual row
  value is at most its formal row value, and the largest is 3.506698068.
- **Damage tests** (`logs/upperB_damage.log`). Each damage raises the exact maximum above R_+, and each is rejected:
  - U of the argmax window −0.1 %;
  - E at its cut ±0.1 %, keeping the symmetries;
  - F ±0.1 %;
  - G ±0.1 %;
  - a zero F is rejected at load.

**Gaps (paper level).**
1. The certificate encoding is undocumented. Paper B §7 says the file holds "the tables as exact rationals, and R_+";
   in fact R_+ is in a separate file (`r9_run2_verify.json`), the candidate's status field says
   NUMERICAL_CANDIDATE_REQUIRES_EXACT_REPLAY, and the lookup conventions are stated nowhere.
2. Paper B does not say whether the padding letters belong to the left and right sides whose grams enter W. By
   Lemma 8 the ratios do not depend on this, so it is harmless, but the definition of W is incomplete as written.
3. Proposition 9 is proved only "in idea". My actual-tape test supports it, and I found no counterexample.

## 4. Paper B, lower bound (Theorem F, Section 6): NOT CHECKED

**Why it could not be checked within the rules.**
- The certificate is a chain of five stages: `1204_coarse_coordinates/data/alt7/*` and `1228_combined_lower/runs/*`.
  It consists of:
  - a radius-6 exact parent stored only as a Python pickle (`exact_parent_physical_certificate.pkl`); opening it
    would run project classes;
  - 968,894 closed-walk "directions" stored as pairs of encoded int64 descriptors, plus float64 coefficients.
- Neither Paper B nor any folder document describes how a descriptor encodes an anchor and a walk, how windows with
  heights are keyed, or the parent's layout. Decoding them would need the producer or checker code, which this brief
  forbids.
- Even with formats in hand, a check from Paper B alone would need the lifting, closed-walk and sparse-refinement
  constructions as precise statements. Paper B gives them as prose in §6.2, not as numbered results.
- The balance equations at radius 9 and 10 are, by Paper B §7, not checked one by one anywhere (0.73 % sampled).
  The bound therefore rests on those constructions, not on a direct check of (S).

**What I did check** (`lowerB/hp_lawB.py`, `logs/lowerB_law.log`). I rebuilt the memory-three reference law from
Paper B §6.1 alone, in exact rationals. The text needs two interpretive choices:
- 4- and 5-paths that leave the strip count as probability 0 in the "smallest probability";
- the orbit is that of the middle heights under reversal and reflection.

Without the first choice the claim "neither change alters shorter paths" fails at the boundary. Under this reading:
- there are 154 five-paths, all positive (minimum 7.5·10⁻¹¹);
- all 3-, 4- and 5-path marginals are consistent on both sides;
- the law is reflection- and reversal-symmetric;
- the memory-three chain is irreducible.

So the hypotheses of Theorem F on the law hold for the law as the paper describes it. I could not compare it with
the project's law, which is stored only in the pickle. A side note: the stored metadata in `combined_cycle_r8.npz`
has `"memory_eta": "1/8"` next to `"memory_parameters": ["3/16","1/8"]`. The name looks wrong, since Paper B has
η = 3/16 and ζ = 1/8.

**What would make bound 4 checkable.**
- A written format for the descriptors and the parent (a non-pickle export).
- Numbered statements of the lifting and closed-walk propositions (1252 Propositions B–D) in Paper B.
- Ideally, an exported final radius-10 flow, or a per-stage list of changed (window, move, value) triples, so that
  (S) can be checked directly in C with exact arithmetic.
