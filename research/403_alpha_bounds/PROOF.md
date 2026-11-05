# 403: capped blocks, a pointwise bound of order 1/n, and the exponent bound 3/2

> **Referee result (2026-10-07, `404_referee_alpha_bounds/REPORT.md`): Theorems 1 to 6 are correct and Theorem 7 is correct as a conditional statement.** What to read differently below:
> 1. Theorem 1's two planarity steps were argued by a picture; the referee wrote them out combinatorially (its section 2). With that repair the bijection is proved.
> 2. Theorems 3 and 5 are proved outright: the existence of μ (input I1) is a theorem (supermultiplicativity plus M_n ≤ 16^n).
> 3. Strengthening: T(A) ≥ 4a − 1, hence **M_n ≤ μ^(n+1)/(8n)** for every n ≥ 1.
> 4. "A block with more than half the points is unique" is false as worded; the correct statement is that a block of a given size 2b+2 ≥ N+1 is unique (blocks with more than N points form a nested chain).
> 5. Theorem 5 is a limsup statement: for every κ < 3/2, liminf n^κ M_n μ^(−n) = 0. It gives α ≥ 3/2 only if a power law exists; it is not a bound at every n (that bound has exponent 1).
> 6. Hypothesis H is open and is about as strong as what it would prove.

All statements in sections 1 to 6 are **proved** here, with the inputs named. Section 7 is a
conditional theorem. Exhaustive small-order checks are listed in section 8. Nothing in this file
needs the existence of a power law unless it says so.

Notation. A meander of order n >= 1 is a pair C = (U, Lo) of noncrossing perfect matchings of the
points 1, ..., 2n on a line, U drawn as arcs in the upper half-plane and Lo in the lower half-plane,
whose union is one closed curve. U(x) and Lo(x) are the partners of x. M_n is the number of meanders
of order n, M_0 = 1, and mu = lim M_n^(1/n). Gap i (0 <= i <= 2n) lies between points i and i+1;
gaps 0 and 2n are the two ends. Write m_n = M_n mu^(-n).

Inputs taken from outside this folder:

* (I1) the limit mu = lim M_n^(1/n) exists (classical; it follows from M_(a+b) >= M_a M_b, which is
  the junction product of 120 Theorem F and is reproved in Lemma 5 below).
* (I2) only for the remarks marked so: sum_(n>=1) m_n <= 1/4 (120 Theorem A; internally refereed).

## 1. Capped blocks

**Definition.** Let C = (U, Lo) be a meander of order N. A *capped block* of C is a pair of points
L < R such that

1. {L, R} is an arc of Lo,
2. R >= L + 3,
3. U maps the interval [L+2, R] into itself.

Its *size* is s = R - L + 1.

**Lemma 1.** Let (L, R) be a capped block of size s of a meander C of order N.

(a) Lo maps [L, R] into itself, s is even, s = 2b + 2 with b >= 1.

(b) If s < 2N then U(L) and U(L+1) lie outside [L, R]. If s = 2N then U(L) = L+1.

(c) The graph induced by C on [L, R] is a single path from L to L+1 that visits every point of
[L, R]; it starts L, R, U(R), ... .

*Proof.* (a) Lo is noncrossing and contains {L, R}, so no arc of Lo leaves [L, R]. A perfect matching
on [L, R] needs s even, and s >= 4. (b) By condition 3 the partners U(L), U(L+1) are not in [L+2, R].
So either U(L) = L+1 or both partners are outside [L, R]. If U(L) = L+1 then [L, R] is closed under U
and Lo, so it carries a whole closed curve; C is one curve, hence [L, R] is everything and s = 2N.
Conversely if s = 2N there is no outside. (c) In the induced graph every point has degree 2 except L
and L+1, which have degree 1 when s < 2N (their upper arcs leave). A graph of maximum degree 2 with
exactly two vertices of degree 1 is one path between them plus cycles, and a cycle would be a second
closed curve of C. If s = 2N, remove the arc {L, L+1} from the single cycle. The path leaves L by its
lower arc, which goes to R. ∎

**Lemma 2 (laminar).** Two distinct capped blocks of one meander occupy intervals that are nested or
disjoint. In particular two distinct capped blocks of equal size are disjoint, and two capped blocks
that each have more than N points are nested.

*Proof.* The intervals are spanned by two distinct arcs of the noncrossing matching Lo, which share no
endpoint. Two intervals with more than N of the 2N points intersect. ∎

## 2. Insertion and extraction

For a meander A of order a >= 1 and a gap i of A, let f(A, i) be the face of the upper half-plane
minus the arcs of U_A whose closure contains gap i, and let d_A(i) be the number of arcs of U_A on the
boundary of f(A, i). Put

    T(A) = sum_(i=0)^(2a) d_A(i),      T_a = sum over meanders A of order a of T(A),      T_0 = 1.

**Extraction.** Let (L, R) be a capped block of size 2b+2 < 2N of C, and a = N - 1 - b >= 1. Put
j = U(R), p = U(L), q = U(L+1). Define

* B on the points L+1, ..., R-1 (relabeled 1, ..., 2b in order): Lo_B is the restriction of Lo; U_B is
  the restriction of U to [L+2, R-1] minus j, together with the arc {L+1, j}.
* A on the points outside [L, R] (relabeled 1, ..., 2a in order): Lo_A is the restriction of Lo; U_A
  is the restriction of U minus the points p, q, together with the arc e = {p, q}.
* i = L - 1, the gap of A where the block was.

**Insertion.** Let A be a meander of order a >= 1, i a gap of A, e = {p, q} an arc of U_A on the
boundary of f(A, i), and B a meander of order b >= 1 with j = U_B(1). Put the 2b+2 new points
L = i+1, ..., R = i+2b+2 into gap i, with B on L+1, ..., R-1. Let Lo_C be Lo_A, Lo_B and the arc {L, R}.
Let U_C be U_A without e, U_B without {1, j}, the arc {j, R}, and the two arcs that join {L, L+1} to
{p, q} without crossing each other.

**Theorem 1 (bijection).** Fix a, b >= 1 and N = a + b + 1. Insertion is a bijection from the
quadruples (A, i, e, B) onto the pairs (C, (L, R)) with C a meander of order N and (L, R) a capped
block of C of size 2b + 2. Its inverse is extraction. Consequently

    sum over meanders C of order N of kappa_b(C) = T_a M_b,

where kappa_b(C) is the number of capped blocks of size 2b+2 of C. Also, the meanders of order N that
have a capped block of size 2N are in bijection with the meanders of order N-1 (so the formula holds
for a = 0 with T_0 = 1).

*Proof.* *Extraction is well defined.* By Lemma 1(a) Lo_B and Lo_A are noncrossing perfect matchings.
By condition 3 and Lemma 1(b), U_B and U_A are perfect matchings. U_B is noncrossing: an arc {x, y} of
U with L+1 < x < j < y < R would cross {j, R}. U_A is noncrossing: the arc {p, L}, the segment
[L, L+1] of the line and the arc {L+1, q} form a simple path from p to q in the closed upper
half-plane that meets no other arc of U; after deleting the points of [L, R] and their arcs, push the
path off the line to get an arc from p to q in the open upper half-plane disjoint from the remaining
arcs. The same picture shows that e is on the boundary of the face of U_A that contains the gap where
[L, R] was. By Lemma 1(c) the closed curve C is the path L, R, j, ..., L+1 through all of [L, R]
followed by a path from q to p through all outside points. Closing the first path minus its first two
points by {L+1, j} gives one closed curve B through L+1, ..., R-1 (here j != L+1 because b >= 1), and
closing the second by {p, q} gives one closed curve A (p != q because a >= 1).

*Insertion is well defined.* Lo_C is noncrossing: every arc of Lo_A over gap i passes over the whole
new interval, and {L, R} passes over B. For U_C: among the two pairings of {L, L+1} with {p, q}
exactly one is noncrossing (if p, q are on the same side of the gap the farther one goes to the
farther of L, L+1; if they are on opposite sides each goes to the nearer one). Because e is on the
boundary of f(A, i), the two new arcs can be drawn inside f(A, i) along e, so they cross no arc of
U_A. All remaining arcs of the block have both ends in [L+2, R] and lie under every arc that passes
over gap i and under any new arc that leaves L or L+1 to the right. Finally {j, R} crosses no arc
{x, y} of U_B with 1 < x < j < y, since that arc would cross {1, j} in B. The curve: A minus e is a
path from p to q through all points of A, B minus {1, j} is a path from j to L+1 through all points
of B, and the new arcs join them through L and R into one closed curve. The pair (L, R) satisfies
conditions 1 to 3.

*The two maps are inverse.* Extraction after insertion returns (A, i, e, B) by construction.
Insertion after extraction returns C because the two upper arcs at L and L+1 of C form a noncrossing
pairing of {L, L+1} with {p, q}, and that pairing is unique.

*Size 2N.* If (1, 2N) is a capped block then U(1) = 2 and the extraction of B alone (as above, with no
A) gives a meander of order N-1; conversely for B of order N-1 the same construction with the arc
{L, L+1} in place of the two joining arcs gives a meander wrap(B) of order N. These are inverse. ∎

**Lemma 3.** For every meander A of order a >= 1, T(A) >= (4a^2 + a + 1)/(a + 1) >= 4a - 3.

*Proof.* U_A has a arcs and a + 1 faces. A bounded face under an arc with k arcs directly below it has
k + 1 arcs on its boundary and contains k + 1 gaps. The unbounded face has c arcs on its boundary (the
exterior arcs) and contains c + 1 gaps (including both ends). So T(A) = sum_f c_f^2 + c_infinity,
where c_f is the number of arcs on face f. Each arc is on two faces, so sum_f c_f = 2a, and by
Cauchy-Schwarz sum_f c_f^2 >= 4a^2/(a+1). Also c_infinity >= 1. ∎

## 3. Theorem 2: insertion inequality with a linear factor

**Theorem 2.** For 1 <= a <= b,

    M_(a+b+1) >= T_a M_b >= ((4a^2 + a + 1)/(a + 1)) M_a M_b >= (4a - 3) M_a M_b.

For all a, b >= 1, T_a M_b <= floor((a+b+1)/(b+1)) M_(a+b+1).

*Proof.* Let N = a + b + 1. Distinct capped blocks of size 2b + 2 of one meander of order N are
disjoint (Lemma 2), so there are at most floor(2N/(2b+2)) of them, and at most one when
2(2b+2) > 2N, that is when a <= b. Apply Theorem 1 and Lemma 3. ∎

Remark (what the factor is). The 2a + 1 gaps of A and the choice of the arc e give T(A) ways to hang
a block made from B on A. The inequality is unconditional because the inserted block is at least as
large as A: a second block of that size cannot fit. This is the device of Madras (J. Stat. Phys. 78
(1995) 681) for lattice trees. Folder 120, section 6, records the same insertion without the size
restriction and notes that its multiplicity cannot be bounded; the restriction a <= b removes the
multiplicity.

## 4. Theorem 3: a pointwise bound of order 1/n

**Theorem 3.** For every n >= 1,

    M_n <= mu^(n+1) / (8n - 2).

*Proof.* Put h(n) = (8n - 2) M_n. Theorem 2 with a = b = n gives
M_(2n+1) >= ((4n^2 + n + 1)/(n + 1)) M_n^2, so

    h(2n+1) = (16n + 6) M_(2n+1) >= ((16n + 6)(4n^2 + n + 1)/(n + 1)) M_n^2 >= (8n - 2)^2 M_n^2 = h(n)^2,

because (16n + 6)(4n^2 + n + 1) - (8n - 2)^2 (n + 1) = 8n^2 + 50n + 2 > 0. Let n_0 = n and
n_(k+1) = 2 n_k + 1, so n_k + 1 = 2^k (n + 1). Then h(n) <= h(n_k)^(2^(-k)) for every k. As
k -> infinity, (8 n_k - 2)^(2^(-k)) -> 1 and M_(n_k)^(2^(-k)) = (M_(n_k)^(1/n_k))^(n_k / 2^k) ->
mu^(n+1) by (I1). ∎

So m_n <= mu/(8n - 2) for every n: an unconditional pointwise bound with exponent 1. If a power law
M_n ~ A mu^n n^(-alpha) holds this gives alpha >= 1, which is weaker than the known alpha > 1; the
content of Theorem 3 is that it holds at every n without regularity assumptions. (The previous
unconditional pointwise bound in this project was of order n^(-1/2), RIGOROUS_RESULTS item 13; the
bound mu^n/(4n) of 120 is conditional on M_(n+1) <= mu M_n.)

Corollary (proved modulo I2). sum_(n>=1) ((4n^2 + n + 1)/(n + 1)) m_n^2 <= mu sum_n m_(2n+1) <= mu/4.

## 5. Theorem 4: a three-piece inequality

**Lemma 4 (wrap).** For a meander B of order r >= 1 let wrap(B) be the meander of order r + 1 of
Theorem 1 (size 2N case): points 1, ..., 2r+2, B on 2, ..., 2r+1, lower arc {1, 2r+2}, upper arcs
{1, 2} and {j, 2r+2} replacing the upper arc {2, j} of B at its first point. B -> wrap(B) is injective
and wrap(B) has the lower arc {1, 2r+2}.

**Lemma 5 (junction product; this is 120 Theorem F).** For meanders B_1 of order p and W of order w,
let B_1 . W be the pair of matchings on 2p + 2w points obtained by placing B_1 before W, keeping all
lower arcs, and replacing the upper arcs {x, 2p} of B_1 and {2p+1, y} of W by {x, y} and {2p, 2p+1}.
Then B_1 . W is a meander of order p + w, and (B_1, W) is determined by B_1 . W and p.

*Proof.* The arc {x, 2p} is exterior in B_1 and {2p+1, y} is exterior in W, so {x, y} passes over
everything between and {2p, 2p+1} over nothing; no crossing appears. Re-pairing one arc of each of
two closed curves gives one closed curve. Given p, re-pair {x, y}, where y = U(x) is the partner of
the unique x <= 2p - 1 with U(x) >= 2p + 2, and {2p, 2p+1} back. ∎

**Theorem 4.** For N >= 5 let R_N be the set of triples (a, p, r) of positive integers with

    a + p + r = N - 2,   2a <= N - 1,   2p <= N,   2r + 2 <= N.

Then

    M_N >= sum over (a, p, r) in R_N of T_a M_p M_r.

*Proof.* For (a, p, r) in R_N, meanders A, B_1, B_2 of orders a, p, r, a gap i of A and an arc e on
f(A, i), let B = B_1 . wrap(B_2), of order b = p + r + 1 = N - 1 - a, and let C be the insertion of
(A, i, e, B), of order N, with capped block (L, R) of size 2b + 2. We show that C determines all the
data, including (a, p, r).

Since 2a <= N - 1 the block has 2N - 2a >= N + 1 points. Capped blocks with more than N points are
nested (Lemma 2). Let (L', R') be a capped block of C strictly inside [L, R]. Then {L', R'} is a lower
arc of C inside [L, R] other than {L, R}, hence a lower arc of B. The lower arcs of B = B_1 . wrap(B_2)
are those of B_1, inside its first 2p points, and those of wrap(B_2), inside its last 2r + 2 points.
So R' - L' + 1 <= max(2p, 2r + 2) <= N. Hence (L, R) is the smallest capped block of C with more than
N points, and C determines (L, R), then a, then (A, i, e, B) by Theorem 1. The last point of B is the
right end of the lower arc {2p+1, 2b} of wrap(B_2), so B determines p, then (B_1, wrap(B_2)) by
Lemma 5, then B_2 by Lemma 4. ∎

## 6. Theorem 5: the exponent is at least 3/2

**Theorem 5.** For every kappa < 3/2,

    liminf_(n -> infinity) n^kappa M_n mu^(-n) = 0.

Equivalently limsup_n log(mu^n / M_n) / log n >= 3/2. In particular, if M_n = mu^n n^(-alpha + o(1))
(for instance if M_n ~ A mu^n n^(-alpha) with A > 0), then alpha >= 3/2.

*Proof.* Suppose not: there are kappa < 3/2, c > 0 and n_0 with M_n >= c mu^n n^(-kappa) for all
n >= n_0. Increasing kappa we may assume 1 <= kappa < 3/2. Let N >= max(16, 4 n_0 + 8). For integers
a, p in [N/4, 3N/8] put r = N - 2 - a - p. Then N/4 - 2 <= r <= N/2 - 2, so a, p, r >= n_0, all three
are at most N/2, and (a, p, r) is in R_N. There are at least (N/8 - 1)^2 >= N^2/256 such pairs (a, p).
By Lemma 3, T_a >= (4a - 3) M_a >= (N/2) M_a. Theorem 4 gives

    M_N >= (N^2/256)(N/2) c^3 mu^(N-2) (N/2)^(-3 kappa) = c_1 mu^N N^(-(3 kappa - 3)),

with c_1 = c^3 8^kappa / (512 mu^2). So the hypothesis holds with kappa replaced by
kappa_1 = 3 kappa - 3 and (c, n_0) replaced by (c_1, max(16, 4n_0 + 8)). Since
3/2 - kappa_1 = 3 (3/2 - kappa), after finitely many repetitions (each one starting from an exponent
in [1, 3/2)) we reach an exponent kappa_j < 1 with M_n >= c_j mu^n n^(-kappa_j) for all large n. This
contradicts Theorem 3. ∎

Remarks.

* What it says. There is no eventual lower bound M_n >= c mu^n n^(-kappa) with kappa < 3/2. It does
  not say that M_n <= C mu^n n^(-3/2) for all n; the bound valid at every n is Theorem 3.
* Existence of a power law is not assumed in the liminf statement and is assumed in the last sentence.
* Why 3/2 is the limit of this method. The inequalities of Theorems 2 and 4 only use insertion of
  blocks and the junction product. A family of objects closed under two such operations can have a
  generating function with a square-root singularity (exponent exactly 3/2), and Theorem 4 is then
  saturated: with m_n of order n^(-3/2) the right-hand side is of order N^2 . N . N^(-9/2) =
  N^(-3/2). The depth-capped meanders of 268 have exponent 3/2. So no argument that uses only these
  operations can prove more; see section 7 for the missing ingredient.

Corollary (proved modulo I2). sum_N mu^(-N) sum_(R_N) T_a M_p M_r <= 1/4, that is
sum over all (a, p, r) in R_(a+p+r+2) of (T_a mu^(-a)) m_p m_r <= mu^2/4.

## 7. Conditional theorem: what would give exponent 2

For N >= b >= 1 let F(b, N) be the number of meanders of order b with no capped block of more than N
points, and F(b, N) = M_b if 2b <= N. Let

    gamma_n = F(n, n) / M_n,

the fraction of meanders of order n with no capped block of more than n points (more than half).

**Theorem 6 (exact identity; proved).** For N >= b >= 2,

    M_b = F(b, N) + sum_(a=0)^(floor(b - (N+1)/2)) T_a F(b - 1 - a, N).

*Proof.* A meander B of order b either has no capped block with more than N points, or has one; in
the second case these blocks are nested (Lemma 2, since N >= b) and there is a smallest, of size
2b' + 2 >= N + 1, with a = b - 1 - b'. By Theorem 1 the pairs (B, smallest block) with given a
correspond to the (A, i, e, B') for a >= 1, or to B = wrap(B') for a = 0, such that no capped block
of B with more than N points lies strictly inside the inserted block. A capped block (L', R')
strictly inside [L, R] has its lower arc in B', and [L'+2, R'] lies in [L+3, R-1]; U_B and U_(B') agree
there except at j, whose partner is R in B and the first point of B' in B', both outside [L'+2, R'].
So the capped blocks of B strictly inside [L, R] are exactly the capped blocks of B' (the whole of B'
included when B' = wrap(B'')). Hence the condition is that B' is counted by F(b', N). ∎

**Hypothesis H.** liminf_n gamma_n > 0.

**Theorem 7 (proved modulo H).** Assume H. Then for every kappa < 2,
liminf_n n^kappa M_n mu^(-n) = 0; in particular a power law M_n = mu^n n^(-alpha + o(1)) forces
alpha >= 2. If in addition M_(n+1) <= mu M_n for all n (hypothesis R; true for n <= 39 by 120), then

    sum_(n >= 1) (4n - 3) M_n mu^(-n) <= mu / gamma_*  < infinity,   gamma_* = inf over n of gamma_n,

that is, the first critical moment is finite (target T1).

*Proof.* Take b = N in Theorem 6 and drop F(N, N) and the term a = 0:
M_N >= sum_(a=1)^(floor((N-1)/2)) T_a F(N-1-a, N). Since N >= N-1-a, F(N-1-a, N) >= F(N-1-a, N-1-a)
= gamma_(N-1-a) M_(N-1-a). Let gamma > 0 and n_1 be such that gamma_n >= gamma for n >= n_1. Then

    M_N >= gamma sum_(a=1)^(min((N-1)/2, N-1-n_1)) T_a M_(N-1-a).                      (*)

If M_n >= c mu^n n^(-kappa) for n >= n_0 with 1 <= kappa < 2, restrict (*) to a in [N/4, N/2 - 1]:
M_N >= gamma (N/4 - 1)(N - 3) c^2 mu^(N-1) (N/2)^(-kappa) N^(-kappa) >= c_1 mu^N N^(2 - 2 kappa) for
large N. The exponent map kappa -> 2 kappa - 2 has the repelling fixed point 2, so finitely many
steps give an exponent below 1, contradicting Theorem 3. Under R, M_(N-1-a) >= mu^(-a-1) M_N, and
(*) with n_1 = 1 gives 1 >= (gamma_*/mu) sum_(a <= (N-1)/2) T_a mu^(-a); let N -> infinity and use
Lemma 3. Under H the infimum gamma_* is positive because gamma_n > 0 for every n: the junction
product of n copies of the meander of order 1 has only lower arcs {2k-1, 2k}, hence no capped block. ∎

Status of H. **Not proved.** Data (section 8, exact values): gamma_2 = gamma_3 = 0.5, then
0.714, 0.702, 0.767, 0.761, 0.794, 0.790, 0.810, 0.809, 0.822 (n = 4, ..., 12), reaching 0.8525 at
n = 24 and 0.8526 at n = 25; the minimum over n <= 25 is 0.5. Heuristic only: if the first moment is
finite one expects gamma_n -> 1/(1 + Theta) with Theta = sum_(a >= 0) T_a mu^(-a-1); the partial sum
to a = 12 is 0.1228, giving about 0.89. Heuristic only: H should fail for tree-like families
(exponent 3/2), where most of the object lies inside nested blocks; proving H means proving that
meanders are not of that kind, and Theorems 1 to 5 cannot do it, because such families satisfy
inequalities of the same form. ATTEMPT.md gives the reformulation of H by rotation and the routes
that were tried.

## 8. Checks (data)

`verify_insertion.py` (Python, exhaustive):

* Theorem 1 for every (a, b) with a + b + 1 <= 8: every insertion is a meander, the inserted pair is a
  capped block, extraction inverts insertion, the set of images equals the set of all (meander, capped
  block of that size) found by direct search, and its size is T_a M_b.
* Lemma 2 for all meanders of order <= 8; at most one block when a <= b; Lemma 3 for a <= 7;
  the size-2N case (wrap).
* Theorem 4 for N <= 11: all constructed meanders are distinct and the inserted block is the
  smallest capped block with more than N points.
* Lemma 5 for orders up to (3, 3).

`census.c` (exhaustive enumeration, orders <= 12) and `analyze.py`:

* T_a = 3, 15, 98, 729, 5918, 51038, 460410, 4299963, 41283770, 405387990, 4055897970, 41225211400
  (a = 1, ..., 12); T_a/(a M_a) = 3, 3.75, 4.08, 4.34, 4.52, 4.65, 4.76, 4.84, 4.91, 4.97, 5.02, 5.06.
* The count identity of Theorem 1, sum_C kappa_b(C) = T_a M_b, for all N <= 12 and all b.
* Theorems 2 (a <= 12, a + b + 1 <= 40; largest ratio T_a M_b / M_(a+b+1) is 3/8 at a = b = 1) and 4
  (5 <= N <= 25; the ratio of the two sides falls from 0.015 to 0.0003) hold with room on exact data.
* Theorem 6: the recursion for F(b, N) against the direct census for b <= 12; gamma_n for n <= 25.

Reproduce (about 7 minutes in total, one thread, under 100 MB):

    clang -O2 -o census census.c
    for k in 1 2 3 4 5 6 7 8 9 10 11; do ./census $k; done > census.txt     # 30 s
    ./census 12 0 104006 > census12_a.txt; ./census 12 104006 208012 > census12_b.txt   # 3.5 min each
    python3 verify_insertion.py 8 11 > verify_insertion.log                 # 6 s
    python3 analyze.py > analyze.log                                        # 1 s
