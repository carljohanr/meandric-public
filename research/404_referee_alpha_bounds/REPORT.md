# 404: referee report on the exponent bounds of 403 (2026-10-07)

Saved by the main session from the referee's final reply.

Folder: `/Users/carljohan/Documents/Code/meandric/research/404_referee_alpha_bounds/` (scripts, logs, `NOTES.md`). Nothing outside it was changed, and no process is left running. I formed my verdicts on the proofs from `PROOF.md` before reading the author's `NOTES.md`, `ATTEMPT.md` and `REPORT.md`.

Labels used: **[code]** verified by my own code; **[proved]** proved by me here; **[read]** read and found correct; **[not checked]**.

## 1. Verdicts

| Claim | Verdict | Confidence | Basis |
|---|---|---|---|
| Lemma 1 | correct (one sentence to add, C2) | 99% | [read], [code] |
| Lemma 2 (laminar) | correct as written | 99% | [read], [code] orders <= 12 |
| Theorem 1 (bijection) | correct with a repair: planarity steps written out in section 2 | 98% | [proved], [code] |
| Lemma 3 | correct as written; holds for every arch system | 99% | [read], [code] orders <= 12 |
| Theorem 2 | correct as written | 98% | [read], [code] |
| Theorem 3 | correct as written; unconditional | 97% | [read], [code] for the algebra and n <= 40 |
| Lemmas 4, 5 | correct as written | 99% | [read], [code] |
| Theorem 4 | correct as written | 97% | [read], [code] N <= 12 |
| Theorem 5 | correct as written in the body; heading overstated (C5) | 96% | [read], bookkeeping [code] |
| Theorem 6 | correct as written (range N >= b >= 2 is needed) | 97% | [read], [code] b <= 12 |
| Theorem 7 | correct as a conditional statement; one constant slip (C6) | 95% | [read] |
| Hypothesis H | open; not a gap in any proof | n/a | data reproduced [code] |
| I1, I2 | used correctly, and I2 only in the two corollaries | 97% | [read] |

## 2. Theorem 1: the two planarity steps, written out [proved]

Conventions: arcs {x, y} and {z, w} cross when x < z < y < w. An arc g = {x, y} is *over* the block (or over gap i) when x < L and y > R.

**Face criterion (F).** The faces of the upper half-plane minus a noncrossing arc system are the nodes of the nesting tree. The face of gap i lies directly under the innermost arc g_0 over the gap (or is the unbounded face), and its boundary arcs are g_0 and the children of g_0. Hence e is on the boundary of f(A, i) if and only if, for every arc g != e: g strictly contains e exactly when g is over gap i.

**Step X (extraction: U_A is noncrossing and e = {p, q} satisfies F).** Let p = U(L), q = U(L+1), both outside [L, R] by Lemma 1(b). Since {p, L} and {q, L+1} do not cross, there are three cases:
- opposite: p < L and q > R;
- left: q < p < L;
- right: R < q < p (mirror of left).

Let g = {x, y}, x < y, be any other arc of U with both ends outside [L, R].

- *g does not cross e, opposite case.* If p < x < q < y, then x < L gives p < x < L < y (crosses {p, L}), and x > R gives L+1 < x < q < y (crosses {L+1, q}). If x < p < y < q, then y < L crosses {p, L}, and y > R gives x < L+1 < y < q (crosses {L+1, q}).
- *g does not cross e, left case.* If q < x < p < y, then y < L crosses {p, L}, and y > R gives q < x < L+1 < y (crosses {q, L+1}). If x < q < y < p, then x < q < y < L+1 crosses {q, L+1}.
- *F, first half.* If g is over the block, then not crossing {p, L} and {q, L+1} forces p and q strictly inside (x, y), so g contains e.
- *F, second half.* If g strictly contains e and is not over the block, the opposite case is impossible. In the left case x < q < p < y < L, which crosses {p, L}.

**Step Y (insertion: the new arcs cross nothing).** Let e = {u, v}, u < v, satisfy F, and let g = {x, y} be another arc of U_A, in the coordinates of C.

- *Opposite case (u < L, v > R), new arcs {u, L} and {L+1, v}.* g crosses {u, L} only if x < u < y < L or u < x < L < y. In the first, y < v, so g crosses e. In the second, g is over the gap, so by F it contains e and x < u: contradiction. {L+1, v} is the mirror.
- *Left case (u < v < L), new arcs {v, L} and {u, L+1}.*
  - g crosses {v, L} with x < v < y < L: g is not over the gap; not crossing e gives x < u, so g strictly contains e, contradicting F.
  - g crosses {v, L} with v < x < L < y: g is over the gap but does not contain e, contradicting F.
  - g crosses {u, L+1} with x < u < y < L: either y < v (crosses e) or g strictly contains e without being over the gap.
  - g crosses {u, L+1} with u < x < L < y: g is over the gap but does not contain e.
- *Block arcs.* An arc with both ends in [L+2, R] cannot interleave with an arc that has no end in [L+2, R]. Inside the block, {j, R} crosses {x, y} only if x < j < y, and that arc would cross {1, j} in B.
- *Uniqueness of the pairing.* For u < v < L, {u, L}, {v, L+1} interleave and {u, L+1}, {v, L} are nested. For u < L, v > R, {u, L}, {L+1, v} are disjoint and {u, L+1}, {L, v} interleave. So exactly one pairing is noncrossing.

The remaining parts of the proof I read and found correct: both pieces are single closed curves by Lemma 1(c); the block conditions hold after insertion; the two maps are inverse; the size-2N case. [read]

Computer confirmation [code]:
- N <= 9: the set of all insertions equals the set of all (meander, capped block) pairs found by direct search, for every (a, b); extraction on every directly found pair returns meanders A, B and an arc e on the right face, and re-insertion returns C.
- N = 10: every insertion is a meander with the claimed block and extraction inverts it (6.8 million quadruples, no failure); the count equals my census count.
- Count identity sum_C kappa_b(C) = T_a M_b, with T_0 = 1: all N <= 12, all b. T_a was computed from the topological definition (no arc separates e from the gap), not from the author's face formula.
- Explicit insertion at N = 11 was not run (too slow in Python); the count identity covers it.

## 3. Theorem 2 and its boundary cases

The proof is correct [read]. For N = a+b+1 the block size is 2b+2, and a <= b is equivalent to 2b+2 >= N+1. Two intervals of at least N+1 points among 2N intersect, so they are nested (Lemma 2), so they cannot have equal size.

| Boundary case | Result | Label |
|---|---|---|
| Size N+1 (N odd, a = b) and N+2 (N even, a = b-1) | at most one block of each size, every meander of order <= 12 | [code] |
| Size exactly N (N even, a = b+1) | two disjoint blocks occur for every even N from 4 to 12 | [code] |
| Any size s | at most floor(2N/s) blocks | [code] |
| Inequality without a <= b | false: (4a^2+a+1)/(a+1) M_a M_1 > M_(a+2) for a = 32..38 | [code], exact M_n |

Example at N = 4: U = {1,6},{2,5},{3,4},{7,8}, Lo = {1,4},{2,3},{5,8},{6,7} is a meander with capped blocks (1,4) and (5,8). So the restriction a <= b is necessary, and the theorem states it correctly.

Lemma 3 is correct [read]: T(A) = sum_f c_f^2 + c_infinity, sum_f c_f = 2a over a+1 faces, then Cauchy-Schwarz. It holds for every arch system; I checked all arch systems of order <= 12, where the minimum of T(A) is exactly 4a-1 [code].

On exact data: T_a M_b <= M_(a+b+1) for all 1 <= a <= b with a+b+1 <= 40 (exact T_a for a <= 12, Lemma 3 beyond); the largest ratio is 3/8 at a = b = 1 [code].

## 4. Theorem 3

Correct [read].
- **Algebra.** (16n+6)(4n^2+n+1) - (8n-2)^2 (n+1) = 8n^2+50n+2 > 0 for all n >= 1 [code]. With the crude factor 4n-3 the step would fail exactly at n = 1, 2; PROOF.md uses the sharp factor (4n^2+n+1)/(n+1) there, so small n are covered.
- **Exact data.** h(2n+1) >= h(n)^2 holds for 1 <= n <= 19; h(3)/h(1)^2 = 44/9 and h(5)/h(2)^2 = 2489/196 [code].
- **Limit.** h(n) <= h(n_k)^(2^(-k)) with n_k + 1 = 2^k (n+1). The polynomial factor tends to 1 and M_(n_k)^(2^(-k)) tends to mu^(n+1). Every n >= 1 starts its own chain, so the bound holds at every n, not only along one chain.
- **Check.** ((8n-2) M_n)^(1/(n+1)) <= 9.494 for n <= 40, below the known rigorous mu >= 11.380, so the bound is true at every tabulated n [code].
- **Input.** Only Theorem 2 with a = b and the existence of mu are used.

Optional strengthening [proved], [code] for the minimum: since the c_f are positive integers summing to 2a over a+1 faces, sum c_f^2 >= 4a-2 and T(A) >= 4a-1 (equality for the rainbow). With h(n) = 8n M_n, h(2n+1) >= 8(2n+1)(4n-1) M_n^2 >= 64 n^2 M_n^2, so M_n <= mu^(n+1)/(8n). The leading constant 8 cannot be improved by this argument.

## 5. Theorems 4 and 5

**Theorem 4: correct [read].**
- The block has 2N-2a >= N+1 points.
- A capped block strictly inside it has its lower arc in B = B_1 . wrap(B_2), hence inside the first 2p or the last 2r+2 points, so it has at most N points.
- Blocks with more than N points form a nested chain, so the inserted block is the smallest of them and C determines it.
- Then a, (A, i, e, B) by Theorem 1; p from the lower partner of the last point of B; (B_1, wrap(B_2)) by Lemma 5; B_2 by Lemma 4.

Computer [code]: for 5 <= N <= 12 (1,839,392 constructions at N = 12) every C is a meander, all are distinct, the inserted block is the smallest capped block with more than N points, and my independent decoder returns all the data. Lemma 4 for r <= 6 and Lemma 5 for orders up to (4, 3) hold as bijections.

Are the side conditions needed? I relaxed each one [code]:

| Relaxed condition | Outcome |
|---|---|
| 2r+2 <= N dropped | collisions from N = 7 (sizes (1,1,3) and (3,1,1) give the same C) |
| 2a <= N-1 weakened to 2a <= N | collisions at N = 8 and 10 |
| 2p <= N dropped | collisions at N = 11 (120 of them) |
| 2p <= N weakened to 2p <= N+2 | no collision for N <= 11; a block inside B_1 has at most 2p-2 points |

So the ranges are essentially sharp for this construction and nothing is overstated. On exact M_n the inequality holds for 5 <= N <= 40, with ratio falling from 0.0115 to 0.00007 [code].

**Theorem 5: correct [read].**
- *Negation.* liminf > 0 gives c > 0 and n_0 with M_n >= c mu^n n^(-kappa) for n >= n_0; raising kappa to at least 1 keeps this true.
- *Ranges.* For N >= max(16, 4n_0+8) and a, p in [N/4, 3N/8]: r lies in [N/4-2, N/2-2], all three are at least n_0 and at most N/2, (a, p, r) is in R_N, there are at least N^2/256 pairs, and 4a-3 >= N/2. I checked these for 16 <= N < 1500 [code].
- *Constant.* c_1 = c^3 8^kappa/(512 mu^2) is right.
- *Direction of T_a.* Only a lower bound on T_a is used, inside a lower bound for M_N. That is the direction needed; no asymptotics of T_a enter.
- *Iteration.* 3/2 - kappa_1 = 3(3/2 - kappa), so finitely many steps give an exponent below 1, which contradicts Theorem 3.
- *Equivalent form.* The limsup form is equivalent to the liminf form taken over all kappa < 3/2.
- *Power law.* "If M_n = mu^n n^(-alpha+o(1)) then alpha >= 3/2" follows: for alpha < kappa < 3/2, n^kappa m_n tends to infinity.
- *Limit of the claim.* The argument does not exclude alpha = 3/2 exactly, and the theorem does not claim it does.

## 6. Theorems 6 and 7, hypothesis H

**Theorem 6: correct [read].** For N >= b, blocks with more than N points are nested, so a smallest one exists. The capped blocks strictly inside the inserted block are exactly the capped blocks of B', because U_B and U_(B') differ only at j, whose partner is outside [L'+2, R'] in both. The summation range is right for both parities of N.

Computer [code]: the identity holds with F(b, N) taken directly from my census for 2 <= b <= 12, b <= N <= 2b+1 (99 cases). It fails for N < b in 61 of 66 cases and for b = 1, so the stated range N >= b >= 2 is needed and correct.

**Theorem 7: correct as a conditional statement [read].**
- (*) follows from Theorem 6 with b = N and F(m, N) >= F(m, m).
- The bootstrap kappa -> 2 kappa - 2 is right up to the constant slip C6.
- Under R, 1 >= (gamma_*/mu) sum T_a mu^(-a) follows, and with Lemma 3 this gives the finite first moment.
- gamma_n > 0 for every n is right: the product of n copies of the order-1 meander has no capped block.

**Is H tested by the data?** The numbers in `analyze.log` are exactly gamma_n = F(n, n)/M_n: by direct census for n <= 12 and by the proved recursion for 13 <= n <= 25. I reproduced all of them (0.5 at n = 2, 3, rising to 0.85248 and 0.85257 at n = 24, 25) [code]. Two cautions:
- H is a liminf statement, and 25 exact values cannot test it.
- H is roughly as strong as its conclusion. Since 1 - gamma_N <= sum_a T_a M_(N-1-a)/M_N, H would follow from a finite first moment plus uniform ratio control. So H is a reformulation of the target, not a reduction to something easier. ATTEMPT.md says this.

## 7. Inputs

- **I1.** M_(a+b) >= M_a M_b is Lemma 5 [read], [code]. With M_n <= C_n^2 <= 16^n, Fekete's lemma gives the limit. PROOF.md omits the upper bound; it is trivial. I1 is therefore a theorem, and Theorems 3 and 5 are proved outright.
- **I2.** sum M_n mu^(-n) <= 1/4 is 120 Theorem A. I read its proof in 120's REPORT.md and the internal referee report REVIEW.md, which found it correct; I did not re-verify it. In 403 it is used only in the two corollaries marked "modulo I2", and there correctly. No theorem depends on it. [read]

## 8. Novelty (limited web search; no priority claim either way)

- **Polynomial factor.** I found no rigorous bound on the sub-exponential factor of meander numbers. The exponent (29+sqrt 145)/12 is a conjecture of Di Francesco, Golinelli and Guitter. The rigorous results I found concern the growth rate (Albert and Paterson, 11.380 <= mu <= 12.901). The higher-genus meander paper of 2023 says the classical asymptotics are still conjectural.
- **Method.** Joining two objects of equal size in order-n ways with a unique way back is Madras's device (J. Stat. Phys. 78 (1995) 681; theta >= (d-1)/d for lattice trees, also 2D polygons). I saw only the abstract.
- **Insertion.** Related operations exist:
  - Lando and Zvonkin's irreducible meander systems. I could not read the PDF. [not checked]
  - Belousov's 2-colored operad on open meanders, which inserts a meander at an intersection point. That paper says the technique does not apply directly to closed meanders; it also notes that iterated snakes have a factor n^(-3/2).
  - The 1-irreducible decomposition M = P/(1-P) in "Meander, folding and arch statistics", the relative of the junction product (from a search summary; I did not read the paper).
  
  I did not find the capped-block insertion, the identity sum kappa_b = T_a M_b, or the inequality M_(a+b+1) >= (4a-3) M_a M_b.

## 9. My computations

All single thread; the longest run was 5 minutes; under 100 MB.

| Check | Range | Result |
|---|---|---|
| M_n, T_n, blocks by size, largest block, multiplicities (`ref_census.c`) | orders 1 to 12 | equal to the author's numbers; T_12 = 41225211400 |
| Laminar property; blocks with more than N points nested | every meander of order <= 12 | 0 violations |
| Count identity of Theorem 1 | N <= 12, all b | pass |
| Explicit bijection, set equality (`ref_bijection.py`) | N <= 9 | pass, 44 cases, 0 bad |
| Explicit insertion, injective plus count | N = 10 | pass |
| Lemma 3 on every arch system | orders <= 12 | pass; minimum 4a-1 |
| Theorem 2, and its failure for a > b | orders <= 40 | pass; failure at a = 32..38, b = 1 |
| Theorem 3 step and bound | n <= 19, n <= 40 | pass |
| Theorem 4 injectivity and decoder (`ref_three.py`) | 5 <= N <= 12 | pass |
| Theorem 4 with relaxed conditions | N <= 11 | collisions as in section 5 |
| Theorem 4 inequality on exact M_n | 5 <= N <= 40 | pass |
| Theorem 6 identity; gamma_n | b <= 12; n <= 25 | pass; equal to the author's values |

Exact M_n for n <= 40 were read from `web/src/meandric/meandricData.ts` and agree with `data/A005315.txt` (29 values) and `data/Jensen_Closed_Meanders.ser` (32 values).

The author's scripts were not audited line by line [not checked]; every number they report that I recomputed agrees with mine.

## 10. Corrections

- **C1 (wording only).** Theorem 1: replace "the same picture shows" and "can be drawn inside f(A, i) along e" by the combinatorial steps X and Y of section 2, and state the face criterion F.
- **C2 (wording only).** Lemma 1(b): add that U restricted to [L+2, R] is a bijection, so no outside point maps into it.
- **C3 (wording only).** I1: mention M_n <= 16^n for Fekete, and label Theorems 3 and 5 "proved" instead of "proved modulo I1".
- **C4 (wording only; 403 REPORT.md, Theorem 2).** "A block with more than half the points is unique" is false as written. Blocks with more than N points form a nested chain and several can occur (wrap(wrap(B)) has two). The correct statement, as in PROOF.md, is that a block of a given size 2b+2 >= N+1 is unique.
- **C5 (wording only).** The heading "Theorem 5: the exponent is at least 3/2" and the folder title should read "upper exponent" or "if the exponent exists". Likewise the remark "why 3/2 is the limit of this method" in section 6 of PROOF.md is a heuristic, while the file header says sections 1 to 6 are all proved; mark it as heuristic there. The reference to "depth-capped meanders of 268" is [not checked].
- **C6 (constant only; does not change a statement).** Theorem 7, proof: the number of integers a in [N/4, N/2-1] is at least N/4-2, not N/4-1 (N = 4k+1 gives k-1). It is absorbed into c_1.
- **C7 (changes a statement; optional strengthening).** Lemma 3 can read T(A) >= 4a-1, and Theorem 3 then reads M_n <= mu^(n+1)/(8n).
- **C8 (wording only).** "Theorem 3 replaces n^(-1/2) by n^(-1)" is true as an order. Numerically the older bound sqrt(pi mu) mu^n/(16 sqrt n) is smaller for n <= 15, and at n = 16 at the estimated mu; Theorem 3 is smaller for n >= 17 (rigorously, using mu <= 12.901).
- **C9 (wording only; 403 REPORT.md).** "H is the single inequality asked for in T4": H implies inequality (*), it is not the same statement.
- **C10 (remark only).** Theorem 4 holds with 2p <= N+2; the other two side conditions cannot be relaxed for this construction.

## 11. What is now proved about the meander exponent, and what is not

Proved, with no unproved input:
- For every n >= 1, M_n <= mu^(n+1)/(8n-2), and in fact mu^(n+1)/(8n). This is a bound of order 1/n at every single n, with no regularity assumption.
- For every kappa < 3/2, n^kappa M_n mu^(-n) comes arbitrarily close to 0 along some sequence of n. Equivalently, limsup log(mu^n/M_n)/log n >= 3/2.
- If M_n = mu^n n^(-alpha+o(1)), then alpha >= 3/2.

Not proved:
- that a power law exists;
- that M_n <= C mu^n n^(-3/2) for all n. The bound valid at every n has exponent 1;
- that alpha > 3/2, or anything about exponent 2. A finite first moment sum n M_n mu^(-n) holds only under the open hypothesis H together with the open ratio hypothesis R;
- any polynomial lower bound on M_n mu^(-n);
- anything near the conjectured value 3.42.

## 12. Not checked

ATTEMPT.md (rotation picture, the "marked wrapping" obstruction, the failed routes); the heuristic remarks on tree-like families and on the limit of gamma_n; the full texts of Madras and of Lando and Zvonkin; 120 Theorem A beyond reading its proof and its internal referee report; explicit insertion at N = 11; the author's scripts line by line.

Sources:
- [Madras 1995, J. Stat. Phys. 78, 681-699 (abstract)](https://springerlink.bibliotecabuap.elogim.com/doi/10.1007/BF02183684)
- [Albert and Paterson, bounds for the growth rate of meander numbers](https://www.cs.otago.ac.nz/research/publications/oucs-2004-08.pdf)
- [Di Francesco, Golinelli, Guitter, Meanders: exact asymptotics](https://arxiv.org/pdf/cond-mat/9910453)
- [Jensen and Guttmann, Critical exponents of plane meanders](https://arxiv.org/pdf/cond-mat/0004321)
- [Di Francesco, Golinelli, Guitter, Meander, folding and arch statistics](https://arxiv.org/pdf/hep-th/9506030)
- [Lando and Zvonkin, Plane and projective meanders (PDF, not readable by my tool)](https://www.labri.fr/perso/zvonkin/Research/proj-mean.pdf)
- [Nica, Free probability aspect of irreducible meandric systems](https://arxiv.org/pdf/1506.01345)
- [Belousov, Prime factorization of meanders](https://arxiv.org/html/2112.10289v5)
- [Belousov, The 2-colored operad on the set of singular meanders](https://arxiv.org/pdf/2206.03453)
- [Delecroix, Goujard, Zograf, Zorich, Higher genus meanders and Masur-Veech volumes](https://arxiv.org/pdf/2304.02567)
