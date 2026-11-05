# SPEC: the certificate for Paper B's lower bound μ ≥ 12.2191

Folder `1421_export_lower_certificate` (prefix IW), written 2026-10-25. This document specifies the certificate
behind the lower bound

  μ ≥ 12.219100857094 = 6109550428547/500000000000

of Paper B (`paper_series/B_tape_bounds`, Section 6, Table 6). It is meant to be read together with Paper B and the
files in `export/`, and nothing else: a reader should be able to write a checker from these two sources alone. The
export files contain no pickle and need no project code to read. Section numbers "Paper B §x" refer to
`paper_series/B_tape_bounds/sections/`. The theory chain is in `1252_referee_lower_theory/PROOF.md` (cited below as
1252 Proposition A to D and Theorem F′).

Here μ is the meander growth rate, μ = lim M_n^{1/n}, with M_n the number of meanders of order n. Throughout,
log is the natural logarithm, and "exact" means rational arithmetic with no rounding.

Contents. 1 The claim and its logic. 2 Tapes, moves, the strip and windows. 3 The reference law. 4 The radius-6
parent. 5 Closed-walk directions. 6 The final flow and its radius structure. 7 What a checker must verify. 8 How the
bound follows. 9 Interpretive choices. 10 What is only sampled. 11 Files and the round trip. Appendix A: reading the
files. Appendix B: the project's native encodings. Appendix C: a worked example.

## 1. The claim and its logic

**The object.** The certificate is a flow J at window radius r = 10 on the tape machine of Paper B §2, confined to the
strip of height K = 12, together with a reference law β of memory m = 3. The flow is

  J = J_P↑ + Σ_d c_d Γ_d,   the sum over 968,894 closed-walk directions d,

where

- J_P is the *parent*: a signed flow at radius 6 with exact rational values (Section 4), and J_P↑ is its lift to radius
  10 (Section 6.1);
- Γ_d is the *field* of the closed walk d (Section 5.2), and c_d is a dyadic rational coefficient (a float64 number
  read exactly).

**The logic.**

1. *Balance (S) and drift (D) hold for J exactly, by construction.* J_P satisfies the balance and drift equations at
   radius 6 exactly (checked on all rows). Lifting preserves them (1252 Proposition A; Paper B §6.2 "Lifting"). Every Γ_d
   satisfies them because d is a protected closed walk (1252 Proposition B; Paper B §6.2 "Closed walks"). The equations
   are linear (1252 Proposition C).
2. *Positivity (P).* Every one of the 45,638,378 flow values J(w, ξ) at radius 10 is strictly positive.
3. *Entropy.* With the total mass M = Σ J > 0, the normalized flow J/M satisfies the hypotheses of Paper B's Theorem F
   (1252 Theorem F′), so μ ≥ e^{2h} with h = (1/M) Σ_{w,ξ} J(w, ξ) log(ν(w)/J(w, ξ)) and ν(w) = Σ_ξ J(w, ξ).
4. *The number.* h ≥ 1.25150018571024385952…, so e^{2h} ≥ 12.219100857094.

**Words used below.** A *flow value* is a number J(w, ξ). The *density* of a flow value is x(w, ξ) = J(w, ξ)/β(w), where
β(w) is the reference probability of the window w. The *density of a window* is ν̄(w) = Σ_ξ x(w, ξ) = ν(w)/β(w).

## 2. Tapes, moves, the strip and windows

### 2.1 Letters, heights, levels

A *word* is a finite string over the two letters `+` and `-`; in sums, `+` counts as +1 and `-` as −1. A *leveled
word* (λ, w) is a word w = w_0 w_1 … w_{k−1} together with an integer λ, its *level*. Its *heights* are the k+1 vertex
values v_0 = λ and v_{i+1} = v_i + w_i. The leveled word *lies in the strip* if every v_i is in {0, 1, …, 11}. Indices
are 0-based everywhere in this document and in the files.

Relation to Paper B (not needed by a checker): for a tape reached from the root, the heights are the gap numbers of
Paper B (2) plus ⌊K/2⌋ = 6, so the root `+[-]+` has heights 5, 6, 5, 6 (Paper B §5.3, "The strip").

### 2.2 Tapes and moves

A *tape* (λ, w, h) is a leveled word with a head index h, 0 ≤ h < |w|. A move is applied only when the head letter
has two neighbors, 0 < h < |w| − 1. The six moves of Paper B Table 2, as string operations:

| move | condition | new word | new head |
|---|---|---|---|
| L | none | unchanged | h − 1 |
| R | none | unchanged | h + 1 |
| S1 | w_h = `+` | w_h replaced by `+-+` | h + 1 |
| S2 | w_h = `-` | w_h replaced by `+--` | h + 1 |
| S3 | w_h = `-` | w_h replaced by `--+` | h + 1 |
| M | w_h = `+` and (w_{h−1}, w_{h+1}) ≠ (`+`, `+`) | the three letters w_{h−1} w_h w_{h+1} replaced by one letter x | h − 1 |

In M, x = w_{h−1} + w_{h+1} + 1 as ±1: (`-`,`-`) gives `-`; (`-`,`+`) and (`+`,`-`) give `+`. The level never changes,
and no move changes the letter sum. On the strip a move is *legal* if the new leveled word lies in the strip. Only S2
and S3 can fail: if v is the height just before the head letter `-`, S2 creates the height v+1 and needs v + 1 ≤ 11,
and S3 creates v − 2 and needs v − 2 ≥ 0.

*Drift vectors* (Paper B (3)), with a and b the numbers of letters left and right of the head:
Δ(L) = (−1, +1), Δ(R) = (+1, −1), Δ(S1) = Δ(S2) = Δ(S3) = (+1, +1), Δ(M) = (−1, −1).

### 2.3 Windows and flows

Fix a radius r and put n = 2r + 1. A *window* at radius r is a leveled word (λ, w) with |w| = n lying in the strip; its
head is the letter w_r. The *legal moves* at a window are: L and R always; S1 if w_r = `+`; S2 if w_r = `-` and
v_r + 1 ≤ 11; S3 if w_r = `-` and v_r − 2 ≥ 0; M if w_r = `+` and (w_{r−1}, w_{r+1}) ≠ (`+`,`+`). Here v_r is the height
before the head letter. A *flow at radius r* assigns a real number to every legal pair (window, move). A *balance row*
at radius r is a leveled word of n + 2 letters lying in the strip. Counts on the strip K = 12
(`roundtrip/counts.json`, by enumeration):

| radius r | letters n | windows | legal pairs (flow values) | balance rows |
|---|---|---|---|---|
| 6 | 13 | 58,450 | 226,079 | 220,110 |
| 8 | 17 | 829,352 | 3,210,357 | 3,125,954 |
| 9 | 19 | 3,125,954 | 12,103,723 | 11,784,514 |
| 10 | 21 | 11,784,514 | 45,638,378 | 44,431,506 |

### 2.4 Centers and completions

The *center at radius r* of a window (λ, w) at radius R ≥ r is the window at radius r formed by the letters
w_{R−r}, …, w_{R+r}, with level λ + w_0 + … + w_{R−r−1}.

A tape (λ_T, T, h_T) *fits* radius r if h_T ≤ r and |T| − 1 − h_T ≤ r. A window (λ, w) at radius r *completes* a tape
that fits radius r if the head of the tape sits on the head of the window, the letters agree and the heights agree:
w_{r−h_T+i} = T_i for 0 ≤ i < |T|, and λ + w_0 + … + w_{r−h_T−1} = λ_T. The windows completing T are T extended by
r − h_T letters on the left and r − (|T| − 1 − h_T) letters on the right, in all ways that stay in the strip.

## 3. The reference law β (memory three)

Paper B §6.1 describes the law in words; this section gives the formulas the certificate actually uses (they are
those of `914_level_memory_gate/edge_model.py`, `1019_overnight_structural/memory_two.py` and `memory_three.py`).
All quantities are exact rationals. `export/law.json` lists the parameters and every table below as exact fractions,
keyed by the list of consecutive heights; a checker should rebuild the tables from the parameters and compare.

**Parameters.** K = 12. Up-probabilities u_0, …, u_11 =
1, 97/100, 99/100, 47/50, 171/200, 69/100, 31/100, 29/200, 3/50, 1/100, 3/100, 0 (so u_{11−i} = 1 − u_i).
κ = 4/5, η = 3/16, ζ = 1/8.

**Layer 1 and 2: the edge chain.** For heights a, b, c in the strip with |a − b| = |b − c| = 1, let P(a, b, c) be the
probability of stepping from b to c after arriving at b from a:
- if u_b ∈ {0, 1}: P(a, b, b+1) = u_b;
- otherwise put o = u_b/(1 − u_b) · κ^{b−a} (odds of stepping up; b − a = ±1) and P(a, b, b+1) = o/(1 + o);
- P(a, b, b−1) = 1 − P(a, b, b+1).

The stationary law of directed edges is ρ(a, b) = C_{min(a,b)} / (2 Σ_{i=0}^{10} C_i), where C_0 = 1 and
C_i = C_{i−1} · P(i−1, i, i+1)/P(i+1, i, i−1) for i = 1, …, 10. It is the unique solution of
ρ(b, c) = Σ_a ρ(a, b) P(a, b, c) with Σ ρ = 1 (22 directed edges). Triples: π(a, b, c) = ρ(a, b) P(a, b, c) (42 entries).
Base quartets: q⁰(a, b, c, d) = π(a, b, c) P(b, c, d) (80 entries).

**Layer 3a: quartets (memory two).** For each quartet (a, b, c, d) in the strip with b and c both in {1, …, 10}, its
*group* is g = min(b, c, 11 − b, 11 − c); let m_g be the smallest q⁰ over all such quartets with group g. Then

  q(a, b, c, d) = q⁰(a, b, c, d) + η · m_g · σ,   σ = (+1 if a < b else −1) · (+1 if d < c else −1),

so σ = +1 exactly when the first step a→b and the last step c→d point in opposite directions. Quartets with b or c in
{0, 11} keep q = q⁰. The memory-two forward kernel is F₂(a, b, c, d) = q(a, b, c, d)/π(a, b, c).

**Layer 3b: quintets (memory three).** Base quintets p⁰(a, b, c, d, e) = q(a, b, c, d) · F₂(b, c, d, e) (154 entries). For
each quintet with b and d both in {1, …, 10}, its group G is the lexicographically smallest of the four triples
(b, c, d), (d, c, b), (11−b, 11−c, 11−d), (11−d, 11−c, 11−b); let m_G be the smallest p⁰ over all such quintets with
group G. Then

  p(a, b, c, d, e) = p⁰(a, b, c, d, e) + ζ · m_G · σ′,   σ′ = (+1 if a < b else −1) · (+1 if e < d else −1).

Quintets with b or d in {0, 11} keep p = p⁰. (The minima m_g and m_G are taken over the unperturbed q⁰ and p⁰.)

**Cylinder probabilities.** For a leveled word (λ, w) in the strip with heights v_0, …, v_k (k = |w|):

  β(λ, w) = ρ(v_0, v_1) if k = 1;  π(v_0, v_1, v_2) if k = 2;  q(v_0, …, v_3) if k = 3;
  β(λ, w) = q(v_0, …, v_3) · Π_{i=4}^{k} p(v_{i−4}, …, v_i)/q(v_{i−4}, …, v_{i−1})  if k ≥ 4.

**Properties (verified exactly; 1252 FINAL, 1253 FINAL item 4, and by the rebuilt tables in `law.json`).**
- All 154 quintets are positive (smallest 7.52·10⁻¹¹), so every leveled word in the strip has β > 0.
- Σ_e p(a, b, c, d, e) = q(a, b, c, d) = Σ_z p(z, a, b, c, d) for every quartet, and Σ_d q(a, b, c, d) = π(a, b, c) =
  Σ_z q(z, a, b, c) for every triple. So β is a stationary Markov chain whose next step depends on the last four
  heights (memory m = 3 steps), and β(w) equals the sum of β over all one-letter extensions of w on either side.
- p and q are invariant under reversal and under reflection h ↦ 11 − h.
- The chain on states of four consecutive heights is irreducible.

**Conditional extensions.** If a leveled word V extends w on the left and the right, β(V)/β(w) is the probability of
the extra letters given w. By the memory-three property it depends on w only through the first three letters of w
(with their heights) for the left part and the last three letters (with their heights) for the right part.

These are the hypotheses of Paper B's Theorem F on the law ("every path of m + 2 heights in the strip has positive
probability and the chain is irreducible", Paper B §5.3) and H1 to H3 of 1252 §0.

## 4. The radius-6 parent

**Definition.** J_P is a flow at radius 6: one exact rational value J_P(w, ξ) = N(w, ξ)/D for each of the 226,079
legal pairs at radius 6, where D is a positive integer of 14,308 bits and each N is a signed integer of 14,252 to
14,298 bits (no common factor with D overall). It is a *signed* flow: 7,708 values are negative and none is zero.
Its densities x_P(w, ξ) = J_P(w, ξ)/β(w) range from about −4.94·10⁶ to 6.18·10⁶.

**Properties the certificate relies on** (Section 7, checks C2 and C3):
- Σ_{w,ξ} N(w, ξ) = D (total mass 1);
- the drift sums are exactly zero: Σ N(w, ξ) Δa(ξ) = Σ N(w, ξ) Δb(ξ) = 0, as integers;
- all 220,110 balance equations (S_6) of Section 7.1 hold exactly.

How the parent was produced (an exact solve over 20,942 coarse coordinates in folder 1204, with the recipe
`1204_coarse_coordinates/data/alt7/exact_parent_exact_recipe.json`) plays no role in the proof; only the numbers do.

**File `export/parent_r6.bin`, format IWPAR001.** All integers are little-endian and unsigned unless stated; there is
no padding.

| offset | size | content |
|---|---|---|
| 0 | 8 | magic, the ASCII bytes `IWPAR001` |
| 8 | 4 | K (= 12) |
| 12 | 4 | r (= 6) |
| 16 | 4 | count (= 226,079), the number of records |
| 20 | 4 | dlen (= 1,789), the number of bytes of D |
| 24 | dlen | D, an unsigned integer, little-endian (least significant byte first) |

Then `count` records follow, each made of a 12-byte head and a magnitude:

| offset in record | size | content |
|---|---|---|
| 0 | 1 | level λ of the window (0 to 11) |
| 1 | 1 | move: 0 = L, 1 = R, 2 = S1, 3 = S2, 4 = S3, 5 = M |
| 2 | 1 | sign: 0 if N ≥ 0, 1 if N < 0 |
| 3 | 1 | reserved, always 0 |
| 4 | 4 | word code of the 13-letter window word: bit 12 (the most significant of 13 bits) is the first letter w_0, bit 0 is the last letter w_12; a 1 bit is `+`, a 0 bit is `-` |
| 8 | 4 | nb, the number of magnitude bytes |
| 12 | nb | abs(N), an unsigned integer, little-endian |

So J_P(w, ξ) = (−1)^sign · abs(N)/D. The file has no bytes after the last record; its size is 406,418,003 bytes.

**Order of the records (canonical order).** Levels 0 to 11 in increasing order; within a level, the window words in
lexicographic order with `+` before `-` (equivalently, decreasing word code); within a window, the legal moves in
the order L, R, S1, S2, S3, M. Every legal pair at radius 6 occurs exactly once and nothing else occurs. A checker
should assert this, since it shows that the parent is defined on all legal pairs.

**File `export/parent_r6_index.tsv`.** One text line per record in the same order, with tab-separated columns
`index level word label sign numerator_bits flow_approx density_approx`, where the word is written with `+` and `-`,
and the last two columns are float64 approximations of J_P and x_P. It is for inspection only; the exact values are in
the binary file.

## 5. Closed-walk directions

### 5.1 A direction

A direction d is a tuple (λ_d, A_d, h_d, ξ_1 … ξ_k, c_d): an *anchor* tape x_0 = (λ_d, A_d, h_d), a sequence of k
moves, and a coefficient c_d. Put x_j = ξ_j(x_{j−1}) for j = 1, …, k (Section 2.2), and write x_j = (λ_d, T_j, h_j)
(the level never changes). The *edges* of d are the k triples (x_{j−1}, ξ_j, x_j). A valid direction satisfies:

1. *closed*: x_k = x_0 (same word and same head index);
2. *legal*: for every j, the head of x_{j−1} is interior (0 < h_{j−1} < |T_{j−1}| − 1), ξ_j applies to it, and the new
   leveled word (λ_d, T_j) lies in the strip;
3. *protected* (memory 3): every T_j has the same first three letters, the same last three letters and the same
   letter sum as A_d (the level is the same by definition);
4. *fits*: every source tape x_{j−1} fits the radius r_s of the stage the direction belongs to (Section 2.4), and
   r_s ≤ 10;
5. *zero drift*: Σ_j Δ(ξ_j) = (0, 0) (automatic for a closed walk, since Δ records the change of the numbers of
   letters left and right of the head; a checker should still assert it).

All 968,894 directions in the export satisfy 1 to 5 (asserted by `scripts/export.py`, and before that by the
project's checkers and by 1253).

### 5.2 The field of a direction

For a window w at any radius R ≥ r_s and a move ξ, the field of d is

  Γ_d(w, ξ) = Σ_{j : ξ_j = ξ and w completes x_{j−1}} β(w) · β(λ_d, A_d) / β(λ_d, T_{j−1}),

in density form Γ_d(w, ξ)/β(w) = Σ_j β(A_d)/β(T_{j−1}) over the same j. If the walk makes the same move from the same
tape twice, both terms count. This is Paper B §6.2 ("the flow value J(w, ξ) grows by z β(T_0) β(w)/β(T)") and the
field of 1252 Proposition B with edge weights f(e) = β(A_d). For edges whose source word has the same length as the
anchor and the same letters, the ratio β(A_d)/β(T_{j−1}) is 1.

The value Γ_d(w, ξ) does not depend on the radius R at which it is evaluated, in the following sense: if R′ > R ≥ r_s,
then Γ_d at radius R′ equals the lift (Section 6.1) of Γ_d at radius R, because a window at radius R′ completes x_{j−1}
exactly when its center at radius R does (x_{j−1} fits R). So every direction may be evaluated directly at radius 10.

### 5.3 The coefficient

c_d is the exact value of a float64 number: c_d = coef_num / 2^coef_pow2, with coef_num a signed integer and
coef_pow2 ≥ 0 an integer. The coefficients were found by numerical optimization; the optimization plays no part in
the proof, and any other coefficients would also give a flow satisfying (S) and (D) (only positivity and the entropy
would change). No coefficient is zero. The largest |c_d| is 3.50·10⁷ (a head cut of stage 0); in later stages
|c_d| ≤ 8.6·10⁴.

### 5.4 Files `export/dirs_<s>_<name>.tsv`

One file per stage s (Section 6.2). Lines starting with `#` are comments. Every other line is one direction, with
ten tab-separated columns:

| column | content |
|---|---|
| index | 0, 1, 2, … within the file, in the order of the project's file |
| level | λ_d |
| anchor | A_d, as a string over `+` and `-` |
| head | h_d, 0-based index into the anchor |
| moves | ξ_1, …, ξ_k separated by commas, from {L, R, S1, S2, S3, M} |
| coef_hex | c_d as Python `float.hex()` text (for example `-0x1.044f38621d012p-9`) |
| coef_num | signed decimal integer |
| coef_pow2 | decimal integer ≥ 0; c_d = coef_num / 2^coef_pow2 exactly, and equals coef_hex |
| native_kind | the project's kind number (only for the round trip; Appendix B) |
| native_key | the project's key (only for the round trip; Appendix B) |

The last two columns are not needed by a checker. The files are plain ASCII with Unix line ends.

### 5.5 The kinds of walk used

Every direction has one of 53 walk shapes (anchor length, head index, move sequence); `export/stages.json` lists
the shapes of each stage with their counts. In words:

- *head cuts* (Paper B "head cycles"): the walk L, R on a tape of 2r_s letters with the head at index r_s (stage 0:
  427,410 at r_s = 8; stage 1: 8,192 at r_s = 9; stage 2: 32,768 at r_s = 10);
- *four-move edit cycles* S, L, M, R and S, R, M, L, where S is the split that applies at the head letter (S1 on `+`;
  on `-`, S2 before L and S3 before R): on anchors of 2r_s − 3 letters with the head in the middle (13 letters in
  stage 0, 15 in stage 1, 17 in stages 2 to 4), and in stages 3 and 4 also on anchors of 18 letters with the head at
  index 9 (S, L, M, R) or 8 (S, R, M, L), the "maximal-context" four-move cycles;
- *six-move cycles* L S1 S3 R M M, L S3 S3 R M M, R S1 S2 L M M, R S2 S2 L M M on anchors of 15 letters, head 7
  (stages 2 to 4);
- *eight-move cycles*: 16 patterns on anchors of 15 letters with head 7 (stages 2 to 4), and in stages 3 and 4 the
  14 further classes of 1211, on anchors of 15 letters (head 8), 17 letters (head 9) or 19 letters (head 10).

Each tape on each walk has at most 20 letters, so every edge source fits radius 10 (and its own stage radius).

## 6. The final flow and its radius structure

### 6.1 Lifting

Let F be a flow at radius r0 and R ≥ r0. The *lift* of F to radius R is

  F↑(W, ξ) = F(center_{r0}(W), ξ) · β(W)/β(center_{r0}(W))

for every window W at radius R and every legal move ξ at W. The legal moves at W are those at its center (legality
reads only the head letter, its two neighbors and the height before the head). In density form the lift changes
nothing: x↑(W, ξ) = x(center_{r0}(W), ξ). Lifting preserves the total mass and the entropy (1252 Proposition D(a)),
and, for r0 ≥ m = 3, the balance and drift equations (1252 Proposition A).

### 6.2 The final flow at radius 10

For every window w at radius 10 and every legal move ξ at w, the certified flow is J(w, ξ) = β(w) · x(w, ξ) with
density

  x(w, ξ) = x_P(center_6(w), ξ) + Σ_d c_d · Σ_{j : ξ_{d,j} = ξ, w completes x_{d,j−1}} β(λ_d, A_d)/β(λ_d, T_{d,j−1}),

where x_P(u, ξ) = J_P(u, ξ)/β(u) is the parent density (Section 4) and the outer sum runs over all 968,894 directions
of the five stage files. Equivalently J = J_P↑ + Σ_d c_d Γ_d, with everything evaluated at radius 10 (Section 5.2).

### 6.3 Stages, the radius-8 completion and the sparse contexts

The flow was built and checked in five stages s = 0, …, 4, with radii r_0 = 8, r_1 = 9 and r_2 = r_3 = r_4 = 10.
Put F_{−1} = J_P (radius 6) and

  F_s = (F_{s−1} lifted to radius r_s) + Σ_{d in stage s} c_d Γ_d   (at radius r_s),

so that F_4 is J (F_4 already lives at radius 10).

- *The radius-8 completion* F_0 is the lift of the parent to radius 8 plus the 542,910 directions of stage 0. It is
  positive on all 3,210,357 pairs at radius 8 and is a certificate by itself (e^{2h} ≥ 12.213777757442).
- *Sparse contexts at radius 9 and 10.* A direction of stage s ≥ 1 changes F only at the pairs (w, ξ) at radius r_s
  where w completes the source tape of one of its edges with move ξ. The union of these pairs over the stage is its
  *changed set*, and the windows that contain a changed pair are its *changed windows* (contexts). Outside its changed
  set, F_s equals the lift of F_{s−1}. The changed sets are parts of all pairs: 621,335 of 12,103,723 at radius 9
  (5 percent), and 6.2 to 10.5 million of 45,638,378 at radius 10 (14 to 23 percent) for each of stages 2 to 4.

Each F_s is positive and is a certificate by itself. Values recorded by the project's checkers (and reproduced by
1253 with independent code, and again here by the round trip of Section 11):

| s | name | r_s | directions | changed pairs | changed windows | e^{2h} ≥ | h ≥ |
|---|---|---|---|---|---|---|---|
| 0 | base (radius-8 completion) | 8 | 542,910 | (all 3,210,357 checked) | (all 829,352 checked) | 12.213777757442 | 1.25128231945632444639964832 |
| 1 | sparse radius 9 | 9 | 24,576 | 621,335 | 372,479 | 12.214596070132 | 1.25131581790761968611843453 |
| 2 | aggregated radius 10 | 10 | 139,264 | 6,165,672 | 3,053,933 | 12.216487956539 | 1.25139325558530462772001740 |
| 3 | increment 1 | 10 | 131,072 | 10,470,289 | 4,596,747 | 12.218343762668 | 1.25146920479395447021016568 |
| 4 | increment 2 | 10 | 131,072 | 10,423,865 | 5,020,763 | **12.219100857094** | 1.25150018571024385952285217 |

The h column is truncated, not rounded; the full 50-digit lower bounds, the mass intervals and the smallest checked
densities are in `export/stages.json` under `recorded`. The smallest density of the final flow is
0.045105124048506763…, and the total mass is M = 0.96713960459460176936685549729213476873473090440078516807…
(exact closed form below, computed by 1253, `1253_referee_lower_computation/runs/mass_exact.json`).

**Exact mass.** Every edge of a direction carries mass c_d β(λ_d, A_d) in total (the windows completing a tape T carry
β(w)/β(T) summing to 1). So M = 1 + Σ_d c_d · k_d · β(λ_d, A_d), with k_d the number of moves of d. This is a cheap
exact check of the export.

## 7. What a checker must verify

### 7.1 The balance equations (S_r), written out

Let F be a flow at radius r ≥ 3 with densities x(w, ξ) = F(w, ξ)/β(w), and put x = 0 for pairs that do not exist.
For a balance row Y = (λ, y_0 … y_{n+1}) (n + 2 letters, in the strip) let u be its inner window, the letters
y_1 … y_n at level λ + y_0. The equation (S_r) at Y is

  x(Y_L, L) + x(Y_R, R) + Σ_{split sources} x(w_S, S) · β(w_S)/β(Y) + Σ_{merge sources} x(w_M, M) · β(V_M)/β(Y) = Σ_ξ x(u, ξ),

where:
- Y_L is the window of the last n letters y_2 … y_{n+1}, at level λ + y_0 + y_1 (the head moved left onto y_{r+1});
- Y_R is the window of the first n letters y_0 … y_{n−1}, at level λ;
- *split sources*: if the three middle letters y_r y_{r+1} y_{r+2} are `+-+`, then w_S is Y with them replaced by `+`
  and S = S1; if they are `+--`, w_S is Y with them replaced by `-` and S = S2; if they are `--+`, replaced by `-`
  and S = S3 (in all cases at level λ, and only if w_S is a window and S is legal there);
- *merge sources*: for each pair (p, q) of letters with p + q + 1 = y_{r+1} (as ±1; one pair (`-`,`-`) if y_{r+1} = `-`,
  two pairs (`-`,`+`), (`+`,`-`) if y_{r+1} = `+`), V_M is the word y_0 … y_r p `+` q y_{r+2} … y_{n+1} of n + 4
  letters at level λ, and w_M is its middle window (letters 2 to n+1 of V_M, at level λ + y_0 + y_1); a term exists if
  V_M lies in the strip.

This is Paper B's (S), Σ_{w,ξ} J(w, ξ) β_{w,ξ}(Y) = ν(Y°) β(Y)/β(Y°), divided by β(Y), with β_{w,ξ}(Y) = β(Y)/β(w)
for L and R (the two fresh letters), 1 for splits and β(V_M)/β(w) for the merge; it is (S_r) of 1252 §0, and it is the
row of the project's matrix (`1019_overnight_structural/memory_three.py`, `build_memory`). (S_r) has one equation per
balance row: 220,110 at radius 6 and 44,431,506 at radius 10 (Section 2.3).

### 7.2 The checks

| check | statement | size | done in the project |
|---|---|---|---|
| C1 law | rebuild ρ, π, q, p from the parameters (Section 3), compare with `law.json`; assert positivity of the 154 quintets, both marginal identities, the two symmetries, irreducibility | 298 entries | exact, full (1252, 1253, 1390 rebuilt it from Paper B) |
| C2 parent rows | (S_6) of Section 7.1 for F = J_P, exactly | 220,110 rows | exact, full (`verify_signed_parent.py`; 1253 rerun; this folder's round trip) |
| C3 parent mass and drift | Σ N = D; Σ N Δa = Σ N Δb = 0 as integers; the parent is defined on exactly the legal pairs at radius 6 | 226,079 values | exact, full |
| C4 directions | conditions 1 to 5 of Section 5.1 for each direction | 968,894 walks | exact, full (project checkers, 1253, `scripts/export.py`) |
| C5 positivity | J(w, ξ) > 0, that is x(w, ξ) > 0, for every legal pair at radius 10 (Section 6.2) | 45,638,378 values | full, organized sparsely (below) |
| C6 entropy | a rigorous lower bound h_lo ≤ h = H/M with H = Σ_{w,ξ} β(w) x(w,ξ) log(ν̄(w)/x(w,ξ)) and M = Σ_w β(w) ν̄(w), all at radius 10; then e^{2 h_lo} ≥ 12.219100857094 | 45,638,378 terms | outward-rounded, full (sparsely) |
| C7 (S_10) directly | (S_10) of Section 7.1 for the final flow J | 44,431,506 rows | **sampled only** (Section 10) |
| C8 (D) directly | Σ_{w,ξ} J(w,ξ) Δa(ξ) = Σ J Δb = 0 at radius 10 | 45,638,378 values | exact by construction; 1253 recomputed it in closed form |

C1 to C6 are what the proof needs (Section 8). C7 and C8 are implied by C2 to C4 and are the direct checks that the
audit asked for; C7 has never been done exhaustively.

**Organizing C5 and C6 sparsely (as the project did).** Take a pair at radius 10 and let s* be the last stage
whose changed set contains its center at radius r_{s*}; its density is that of F_{s*} at that center (or of F_0 at its
radius-8 center if there is no such stage). So C5 is equivalent to: F_0 > 0 on all 3,210,357 pairs at radius 8, and, for
s = 1, …, 4, F_s > 0 on every legal pair of every changed window of stage s, evaluated with the *cumulative* value
(parent plus all stages up to s). This equivalence needs the changed sets to be complete: every pair where some edge of
the stage has its source completed must be in the set (1252 Proposition D(c)). A checker that evaluates x(w, ξ) by
enumerating, for every direction and edge, all completions of the source tape builds complete changed sets by
construction. Similarly H and M change only on changed windows: H(F + δ) − H(F) = Σ_W [ν′ log ν′ − ν log ν] −
Σ [F′ log F′ − F log F] over the changed windows W and changed pairs (in flow units), and M changes by Σ δ.

**Rounding.** The parent densities reach ±6·10⁶ and stage-0 coefficients 3.5·10⁷, while the final densities are
of order 1 (the smallest is 0.045): there is heavy cancellation, and float64 evaluation loses about eight digits. Either
evaluate x exactly (the parent is a rational of about 14,300 bits; directions add dyadics times ratios of β) or use
intervals of at least 128 bits. The project used 160-bit MPFR intervals with directed rounding: enclosures
β ∈ [β₋, β₊], x ∈ [x₋, x₊], ν̄ ∈ [ν̄₋, ν̄₊] give the safe entropy term β₋ x₋ (−log min(1, x₊/ν̄₋)) per pair and the safe
mass term β₊ ν̄₊ per window (Paper B §6.3, "Rounding outward"); then h_lo = (Σ entropy terms)/(Σ mass terms), and
e^{2 h_lo} is rounded down and truncated to 12 decimals.

**Numbers a new checker should reproduce** (from the project's certificates; all reproduced by 1253 with independent
code to within 3.2·10⁻⁴² in h):
- C2: all 220,110 rows have residual exactly 0; 7,708 negative parent values.
- the five stage values of the table in Section 6.3, including h ≥ 1.2515001857102438595228521727409128195625777562276
  at the end;
- M ∈ [0.96713960459460176936685549729213476873473072658955, 0.96713960459460176936685549729213476873473108232162]
  (the project's interval), M = 0.9671396045946017693668554972921347687347309044007851680766… exactly (1253);
- the smallest final density 0.0451051240485067… (a lower bound 0.045105124048506763523197464102908600022095425429597);
- stage 0: 829,352 windows, 3,210,357 pairs, smallest density 0.052478662419946645789…

## 8. How the bound follows

**Theorem F** (Paper B §5.3; proof outline Paper B §5.4; full proof for strips and memory m in 1252 §2 as Theorem F′,
building on `805_entropy_certificates/PROOF.md` §1). Let β be a reference law of memory m on the strip with every path
of m + 2 heights positive and the chain irreducible. Let r ≥ m and let J be a flow at radius r, with total 1,
satisfying (S), (D) and (P): J(w, L) > 0 for every window w. Then μ ≥ e^{2h}, h = Σ_{w,ξ} J(w,ξ) log(ν(w)/J(w,ξ)).

The certificate meets the hypotheses as follows.

1. *The law* (C1): m = 3; the 154 paths of 5 heights are positive; the chain is irreducible (Section 3).
2. *Balance and drift of the parent* (C2, C3): (S_6) and (D) hold exactly for J_P, with Σ J_P = 1.
3. *Lifting* (1252 Proposition A, r0 = 6 ≥ m): a flow, signed or not, that satisfies (S_{r0}) and (D) satisfies
   (S_R) and (D) after the lift to any R ≥ r0, here R = 10. So J_P↑ satisfies (S_10) and (D).
4. *Closed walks* (1252 Proposition B, r = 10 ≥ m): for a protected closed walk whose edge sources fit radius r, the
   field Γ_d of Section 5.2 satisfies (S_r) and (D), for any real coefficient. The proof pushes the measure
   β(· | x_{j−1}) forward along each edge; protection makes the far tape's conditional law the same at every tape of
   the walk, and closedness makes inflow equal outflow at every tape. (C4 checks the hypotheses.)
5. *Signed parent plus fields* (1252 Proposition C): J = J_P↑ + Σ_d c_d Γ_d satisfies (S_10) and (D) by linearity,
   whatever the c_d, and has mass M = Σ J = 1 + Σ_d c_d k_d β(A_d). If J > 0 on every legal pair (C5), then M > 0 and
   J/M is a flow to which Theorem F applies, with h(J/M) = H/M (the entropy is homogeneous of degree one).
6. *Sparse refinement* (1252 Proposition D): lifting preserves mass and entropy; a change supported on a set of pairs
   changes H and M only there; positivity outside the changed sets is inherited from F_0, provided the changed sets are
   complete; chains of stages add cumulatively. This is what makes C5 and C6 cheap; it is not needed by a checker that
   evaluates all 45,638,378 values directly.
7. *The number* (C6): h(J/M) ≥ 1.2515001857102438595…, so μ ≥ e^{2h} ≥ 12.219100857094 = 6109550428547/500000000000.

Steps 3 to 6 are the constructions that Paper B §6.2 gives in prose ("Lifting", "Closed walks", "A signed parent",
"Sparse refinement"). Their written proofs are 1252 Propositions A to D, written by a referee inside the project who
re-derived the chain independently (status "proved"; confidence 94 percent in `ledger/STATE.md`; Paper B Appendix A:
"proved; re-derived independently (1252)"). Theorem F itself rests on
Theorem 3 of Paper B (the tape machine bijection) and on Lemma 4 of Paper B (linear completion), both proved in the
project and read independently there. Nothing has been reviewed outside the project.

**Why r ≥ m matters.** In Proposition A the target window of a split or a merge must keep the first m and last m
letters of the source window; for a split this needs r ≥ m. With r = m − 1 the balance equations at radius r do not
imply those at larger radius (1252 §1, sharpness data). Here r0 = 6 and r = 10, both ≥ 3.

## 9. Interpretive choices

Each item names a point where Paper B's prose, or the project's files, admits more than one reading, and states the
reading the certificate uses. Items 1 and 2 are the two choices that the audit (1390 `AUDIT.md` §4) had to make when
it rebuilt the law from Paper B alone; items 3 to 6 are its other §4 gaps (formats and constructions). In every case the
reading below is the one the project's checkers implement, confirmed by the round trip (Section 11) or by the
conformance test `scripts/spec_conformance.py`.

1. *Paths that leave the strip, in "the smallest probability".* The certificate perturbs only quartets (a, b, c, d)
   with b and c in {1, …, 10} and quintets (a, b, c, d, e) with b and d in {1, …, 10}; the minimum is taken over such
   paths of the same group. This equals the audit's reading ("off-strip paths count as probability 0"): if b or c (b or d)
   is 0 or 11, the orbit contains a path with an end height −1 or 12, so the minimum is 0 and nothing changes; if they
   are interior, every path of the orbit lies in the strip.
2. *The group of a path.* It is the orbit of the middle heights under reversal and reflection: (b, c) for quartets and
   (b, c, d) for quintets. For quartets the project computes it as min(b, c, 11−b, 11−c), which labels the same orbits
   (for c = b ± 1 both determine the unordered pair {b, c} up to reflection). Also: the minimum is over the unperturbed
   values q⁰ and p⁰, and the quintet layer is applied to the already perturbed quartets (Section 3).
3. *How a direction encodes an anchor and a walk.* The export writes each walk explicitly (Section 5.4). The
   project's packed int64 encoding is in Appendix B, for readers of the original files.
4. *How windows with heights are keyed.* A window is (level, word) with the level at the left end; the project's
   integer key is (level << n) | code (Appendix B).
5. *The parent's layout.* Section 4 (export) and Appendix B (the original pickle: a dict of plain Python integers,
   fractions and one numpy array; no project classes).
6. *The constructions as statements.* Lifting (Section 6.1), the field of a closed walk (Section 5.2), the final flow
   (Section 6.2), stages and changed sets (Section 6.3); their proofs are 1252 Propositions A to D (Section 8).
7. *Which letter is the head of an anchor.* Any interior letter; most anchors have it in the middle, but the 18- and
   19-letter "maximal-context" anchors have it off center (head 8, 9 or 10). Only the fit condition matters.
8. *Head cycles.* Paper B §6.2 calls them R L; the certificate writes each as L, R from head index r_s on a tape of 2r_s
   letters. The two are the same closed walk started at its other tape; the field is the same.
9. *Normalization.* Paper B's flows have total 1. The certificate's J has total M ≈ 0.967 (directions add or remove
   mass), and h is H/M. Theorem F is applied to J/M.
10. *(P).* Theorem F asks only J(w, L) > 0. The certificate has J > 0 on every legal pair, which C5 checks.
11. *Fresh letters.* In (S) the letters brought into view are drawn from β conditioned on the window; by the
    memory-three property this is the same as conditioning on its three end letters (Section 3), as Paper B states.
12. *Coefficients.* Binary floating-point numbers read as the exact dyadic rationals they represent (Paper B §6.2).
13. *Signs.* The parent is signed; only the final flow must be positive (1252 Proposition C).
14. *Radius of a direction.* Directions of a stage at radius r_s are evaluated directly at radius 10; this is the same
    flow (Section 5.2).
15. *Consistency checks of the project's checkers that are not hypotheses.* The "parity" interval (mass on windows
    whose head starts at even height minus mass at odd height, for the change made by a stage), the drift intervals
    of about 10⁻⁴⁷ around 0 (the drift is exactly 0 by construction), and the "physical proximity" gate that compares
    the exact increments with the producer's saved floats. None is part of the proof.
16. *A metadata naming slip.* `exact_parent.json` and the npz metadata contain `"memory_eta": "1/8"` next to
    `"memory_parameters": ["3/16", "1/8"]`. The law uses η = 3/16 (quartets) and ζ = 1/8 (quintets), as Paper B says:
    the parent stores parameters = (3/16, 1/8), and the checkers call MemoryThreeLaw(K, up, κ, η = 3/16, ζ = 1/8). The
    field `memory_eta` names the parameter of the last layer and is not read by any checker.
17. *A chain-of-custody gap (1253 §2).* The project's radius-9 file names as its source the optimizer's float file
    `cycle_r8.npz`, not the accepted stage-0 file `combined_cycle_r8.npz` (same 542,910 walks, different coefficient
    representation). The checkers never compare the two; soundness is unaffected because positivity and entropy are
    recomputed over the accepted stage 0. In the round trip of Section 11 the rebuilt radius-9 file points at the
    rebuilt stage 0.
18. *Strip of the completion.* The counted closed walks of Theorem F's proof leave the strip outside the middle part;
    only the stationary process lives on the strip (1252 Theorem F′, remark (i)). Nothing in the certificate depends on
    it.

## 10. What is only sampled, or not checked anywhere

**Only sampled.**
1. *The balance equations of the final flow at radius 9 and 10 (check C7).* They hold by construction (C2, C3, C4 and
   1252 Propositions A to C), and the proof does not need them checked one by one. But no one has checked them
   directly in full. 1253 checked, exactly and with zero residual: 324,436 radius-10 rows of the final flow (0.73
   percent of 44,431,506), 107,568 radius-9 rows (0.9 percent of 11,784,514) and 371,157 radius-8 rows (12 percent of
   3,125,954), including complete sub-families (all rows with fixed central letters) and rows whose sources mix changed
   and unchanged values. A defect confined to fewer than about 10⁻⁵ of the rows could have escaped these samples (1253
   §4). An exhaustive exact check of C7 would close this.
2. *The lifting identity as an identity of numbers* (Proposition A applied to the stored parent): 1253 checked
   Σ over the 16 radius-10 extensions of β·x = β·x at radius 8 on 3,000 flows. It is proved, so this is a test only.
3. *Agreement of the project's gather with an independent scatter* at radius 9 and 10: on 588,147 random pairs
   (1253). At radius 8 it was compared on all 3,210,357 pairs, and the changed sets of all four sparse stages were
   compared in full (equal).
4. *1253's own row operator and tape machine* were validated exactly against the project's matrix only at K = 6,
   r = 4, and by 10⁵ sampled rows at K = 12.
5. *This folder's conformance test* (`roundtrip/spec_conformance.json`): the field formula of Section 5.2, evaluated in
   float64 from the export alone, against the project's own evaluation, on 20,000 sampled pairs per stage, plus 15,000
   to 19,000 sampled pairs per sparse stage outside the changed set (all zero). Largest relative difference
   2.3·10⁻¹²; it tests this document, not the bound.

**Checked in full** (exactly or with outward rounding): the law (C1); all parent rows, the parent mass and drifts
(C2, C3); every direction (C4), including closure, protection, fit and zero drift; positivity of the final flow on
every changed pair of every stage with cumulative values and on all radius-8 pairs (C5, with complete changed sets:
the producer's saved sets equal the checker's enumeration in every part, and equal 1253's enumeration); the entropy
and mass at every stage (C6); the exact mass and drift in closed form over all directions (1253).

**Not done anywhere.** No proof is formalized in a proof assistant. No proof or program has been reviewed outside the
project (the audit 1390 could not check this bound). The direct exhaustive check C7 has not been run. Positivity on the
45,638,378 radius-10 pairs has been checked only in the sparse organization of Section 7.2, never by a direct pass over
all pairs.

## 11. Files and the round trip

### 11.1 The export (`export/`, 507,606,459 bytes)

| file | bytes | SHA-256 | content |
|---|---|---|---|
| `law.json` | 26,551 | `aaa797ef6fe131790ff57e58167b68b0dddf224dda8e418f371f8b5bf7c9aae5` | law parameters and exact tables (Section 3) |
| `parent_r6.bin` | 406,418,003 | `5afa253a51f8161a81c56aa302dfe1aa60a24fe439588dc081ca30612951832a` | the parent, IWPAR001 (Section 4) |
| `parent_r6_index.tsv` | 17,028,785 | `84e16150171a13d7ae7207826f55a284318d660f95b9673c253b0b879876e247` | parent index, floats only |
| `dirs_0_base_r8.tsv` | 44,747,227 | `728f831664fdfc0e4d39eb4b1d3112650082602434804368320d0c933bf2efd3` | stage 0, 542,910 walks |
| `dirs_1_sparse_r9.tsv` | 2,096,812 | `8ed6c12277587c61c321c74442131bb9b74e1800a686e49aff72705a14616c75` | stage 1, 24,576 walks |
| `dirs_2_aggregated_r10.tsv` | 12,555,377 | `21117dfaa886e4d01313cbe080cf1541a92d37509b4574a5cb6eba61e46cee23` | stage 2, 139,264 walks |
| `dirs_3_increment1_r10.tsv` | 12,409,341 | `038275b4190e3e7b7a34be2eaa1d03602eab29842624e42107438f456a03067a` | stage 3, 131,072 walks |
| `dirs_4_increment2_r10.tsv` | 12,305,087 | `d09f8c3657133b2345dc9ad4e1019a29d43ca79bae63647562795419be91a601` | stage 4, 131,072 walks |
| `stages.json` | 18,668 | `3392ebc058130c35d8a6320cddd85d5b9bb31f775f8c519cd6984cfe2ec66309` | stage list, shapes, recorded values, native hashes |

`SHA256SUMS` in the folder lists these and every script and result file. The project files they were exported from,
with their SHA-256, are named in `stages.json` (`native_pickle`, `native_npz`): the parent pickle `f77de744…2fcee`;
stages `325ed6f4…384d`, `54ef75e0…d484`, `366e9aa2…fe7c`, `45500f94…a64a`, `75ae50c3…e019`, the same hashes the project's
certificates record.

### 11.2 The round trip

1. `scripts/export.py` wrote the export from the project's files, replaying and asserting every walk (closed,
   legal, protected, fitting, drift free) and asserting the canonical order of the parent.
2. `scripts/reimport.py` read only the export and rebuilt, by the encoding of Appendix B, the parent (keys,
   numerators, denominator, law parameters) and the descriptor and coefficient arrays of all five stages. All are
   bit-identical to the project's pickle and npz arrays (`roundtrip/reimport_report.json`). It then wrote npz files
   for the project's checkers under `roundtrip/rebuilt/`: arrays from the export; the auxiliary members (metadata, the
   producer's saved changed set and float increments) copied from the original files, with only the metadata fields
   `source` and `source_sha256` pointed at the rebuilt previous stage. The rebuilt stage-0 file is byte-identical to
   the project's `combined_cycle_r8.npz`.
3. `scripts/roundtrip_chain.sh` ran the project's checkers, unchanged (same files, same SHA-256), through
   `scripts/run_checker.py`, whose only change is that `pickle.load` of a `.bin` file returns the parent read from
   `parent_r6.bin`; any other unpickling raises an error. No pickle was opened in any of the eight runs. The parent
   file named in every certificate is therefore the export file itself (SHA-256 `5afa253a…832a`).

| step | checker (unchanged) | result | fields equal to the project's certificate |
|---|---|---|---|
| parent rows | `1019/verify_signed_parent.py` | PASS, all 220,110 rows exact, both drifts 0, 7,708 negative | all 13 non-file fields |
| stage 0 | `1019/signed_memory_cycle_replay.py` | 12.213777757442 | all 21 |
| stage 1 | `1019/verify_sparse_cycle_increment.py` | 12.214596070132 | all 19 |
| stage 2 | `1019/verify_sparse_mixed_paths.py` | 12.216487956539 | all 19 |
| stage 3 | `1228/verify_tmpl.py`, 5 parts | 12.218343762668 | all 19 |
| stage 4 | `1228/verify_tmpl_split.py`, 12 parts in three runs | **12.219100857094** | all 19 |

"Non-file fields" are all fields except wall time and the SHA-256 of input files (which name the rebuilt files);
`roundtrip/compare_certificates.json` lists them. In particular every entropy bound, mass interval, drift and parity
interval, smallest density, count of changed pairs, windows and checked pairs, and the checker's own SHA-256 are
identical. The intermediate state after parts 0 to 3 of stage 4 equals the producer's `inc2_state_0.json` in every
accumulator. The final round-trip certificate `roundtrip/inc2_certificate.json` has SHA-256
`1a15cc63cbdd022c02807d4e57c0042dc29ef79d4fb6a1a0eb701bedd28ed56f`; it states h ≥
1.2515001857102438595228521727409128195625777562276 and μ ≥ 6109550428547/500000000000.

What the round trip shows: the export contains exactly the data the project's checkers certified, and the reading of
Appendix B is the checkers' reading. It is not an independent check of the bound: it reruns the project's own
checkers. Run times: 48 s, 19 s, 21 s, 168 s, 571 s, 3 × 400 s; peak memory 4.5 GB (stage 2).

## Appendix A. Reading the files (a format reader, not a checker)

Reading `parent_r6.bin` in Python (the same steps work in C with a big-integer library):

```
f = open('parent_r6.bin', 'rb'); assert f.read(8) == b'IWPAR001'
K, r, count, dlen = struct.unpack('<IIII', f.read(16)); D = int.from_bytes(f.read(dlen), 'little')
for i in range(count):
    level, move, sign, zero, code, nb = struct.unpack('<BBBBII', f.read(12))
    N = int.from_bytes(f.read(nb), 'little') * (-1 if sign else 1)
    word = ''.join('+' if code >> (2*r - j) & 1 else '-' for j in range(2*r + 1))
    # flow value J_P(level, word, ['L','R','S1','S2','S3','M'][move]) = N / D
assert f.read(1) == b''
```

A line of `dirs_*.tsv` splits on tabs into the ten columns of Section 5.4; the coefficient is
`Fraction(int(coef_num), 2**int(coef_pow2))`, which equals `float.fromhex(coef_hex)` exactly. `law.json` is JSON:
`up` is a list of 12 fraction strings indexed by height; `rho`, `pi`, `quartet` and `quintet` are lists of entries
`[h_0, …, h_k, "num/den"]` (k = 1, 2, 3, 4), keyed by consecutive heights. `stages.json` gives, for each stage, its
file, radius, number of directions, walk shapes with counts, the native file and its SHA-256, and the values recorded
by the project's checker. `scripts/parent_bin.py` is a complete reader of the parent and the law.

## Appendix B. The project's native encodings (for readers of the original files)

Words are coded as integers with the first letter in the most significant bit, `+` = 1: code(w) = Σ_i [w_i = +]·2^{|w|−1−i}.
Moves are indexed L = 0, R = 1, S1 = 2, S2 = 3, S3 = 4, M = 5.

- *Window key* at radius r (n = 2r + 1): (λ << n) | code(w). *Flow key*: (window key << 3) | move index.
- *Parent* `1204_coarse_coordinates/data/alt7/exact_parent_physical_certificate.pkl` (SHA-256 `f77de744…2fcee`): a
  pickled dict with `K` = 12, `r` = 6, `up` (dict height → Fraction), `kappa` (Fraction), `parameters` = (η, ζ) as
  Fractions, `keys` (numpy uint32 array of flow keys in the canonical order of Section 4), `denominator` (int D) and
  `numerators` (list of int N). The only globals it references are `fractions.Fraction` and numpy's array
  reconstruction.
- *Stage files* (npz = zip of npy arrays): `descriptors`, an int64 array of shape (count, 2) with rows (kind, key);
  `coefficients`, float64; `metadata`, a JSON string; for stages 1 to 4 also `changed_keys` (sorted int64 flow keys at
  radius r_s: the producer's changed set) and `density_increment` (float64, the producer's float increments there),
  which the checkers compare with their own enumeration and use as a sanity gate only.
- *Descriptor (kind, key) at stage radius r_s* (n = 2r_s + 1):
  - kind 0 (head cut): key = (λ << (n−1)) | code(T), |T| = n − 1, head r_s, walk L, R;
  - kind 1 and 2 (four-move cycle): key = (λ << m) | code(A), m = n − 4, head m//2 (integer division); kind 1 walk
    S, L, M, R and kind 2 walk S, R, M, L, with S = S1 if A_head = `+`, else S2 (kind 1) or S3 (kind 2);
  - stage 2 only (`family` = `mixed_literal_cycles`), kinds 3 to 22: key = (λ << m) | code(A), m = n − 6, head m//2,
    walk number kind − 3 in the list: L S1 S3 R M M; L S3 S3 R M M; R S1 S2 L M M; R S2 S2 L M M; then the 16
    eight-move patterns L L S1 R R R M L; L L S3 R R R M L; L S1 R R R M L L; L S1 S2 R R M M L; L S1 S3 R R M M L;
    L S3 R R R M L L; L S3 S2 R R M M L; L S3 S3 R R M M L; R R S1 L L L M R; R R S2 L L L M R; R S1 L L L M R R;
    R S1 S2 L L M M R; R S1 S3 L L M M R; R S2 L L L M R R; R S2 S2 L L M M R; R S2 S3 L L M M R;
  - stages 3 and 4 (`family` = `template_trial_cycles`): the metadata list `templates` maps each kind (56 to 73 and
    78 to 101) to an anchor length m, a head h0 and a move list; key = (λ << m) | code(A). All 42 templates have
    `dtype` = `cycle`.
- *Field evaluation in the project's checkers*: kinds 0 to 2 by the inverse map `terms()` of
  `1019_overnight_structural/signed_memory_cycle_replay.py`, the other kinds by a forward scatter over all completions
  of each edge source. Both compute Section 5.2.

## Appendix C. A worked example

Direction 427534 of stage 0 (`dirs_0_base_r8.tsv`): level 0, anchor A = `++++++-++++++` (13 letters, heights
0, 1, 2, 3, 4, 5, 6, 5, 6, 7, 8, 9, 10, 11), head 6 (the letter `-`), moves S2, L, M, R, coefficient
c = `0x1.128915deb359ap-3` = 2414839707572941/2^54 ≈ 0.134050532206328.

- S2 from (A, 6): the word becomes T = `+++++++--++++++` (15 letters, heights 0, 1, …, 7, 6, 5, 6, …, 11), head 7.
- L from (T, 7): head 6.
- M from (T, 6): the letters `+`, `+`, `-` at indices 5, 6, 7 become one `+`; the word is A again, head 5.
- R from (A, 5): head 6, back to the anchor. The level, the first three letters `+++`, the last three letters `+++` and
  the letter sum 11 are the same at every tape; the drift is (+1,+1) + (−1,+1) + (−1,−1) + (+1,−1) = (0, 0).

β(A) ≈ 7.03703736005282·10⁻¹⁵ and β(T) ≈ 1.34736799200025·10⁻¹⁵, so β(A)/β(T) ≈ 5.22280282879949. The density
increments are: c ≈ 0.134051 for (A, 6, S2) and (A, 5, R), and c·β(A)/β(T) ≈ 0.700119 for (T, 7, L) and (T, 6, M).
At radius 8 the windows completing (A, 6) are A with two letters added on each side and the head on index 8: the
prefixes `+-` (level 0) and `--` (level 2), and the suffixes `-+` and `--`; for example the window `+-++++++-++++++--`
at level 0, whose S2 density gains 0.134051 from this direction. At radius 10 the same tape is completed by windows
with four letters added on each side.
