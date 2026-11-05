# 637: proofs (tape form of the label machine, the matrix-product algebra, the no-go theorem, the D = 1 bound)

Status words as in `WAVE_620_COMMON.md`. Everything here rests on Theorem M of `630_referee_band_growth/REPAIRS.md`
section 2 (refereed): for a band Lo ≤ g ≤ Hi with Lo ≤ −1 and Hi ≥ 1, the number M_m of meanders of order m whose
winding function stays in the band equals the number of walks t → t of length 2m − 2 of the label machine. The
machine is quoted in section 0. Nothing else from the project is used.

## 0. The label machine (quoted)

State (h, L, U): h is the winding at the last point; L is the list of labels of the open lower arcs, bottom to top,
with labels in [Lo, Hi − 1]; U the same for the open upper arcs, with labels in [Lo + 1, Hi]. The lower top is
*unmarked* if it equals h and *marked* if it equals h − 2; the upper top is unmarked if it equals h and marked if it
equals h + 2; an empty stack counts as marked. A state with both tops marked is dead and is never entered.
One move: choose d = ±1 with h′ = h + d in the band; d = +1 needs the lower top unmarked, d = −1 the upper top
unmarked. Then on each side either push h′ (if it is an allowed label of that side) or pop the top; a lower pop is
allowed iff the lower top equals h′ − 1, an upper pop iff the upper top equals h′ + 1. Start t = (0, (0), (0)).

## 1. The tape form [proved]

**Definition.** A *tape* is a finite word w = w_1 … w_N over {+, −} (steps +1 and −1) together with a *head*
position i ∈ {1, …, N}; I write … x [s] y … with s = w_i. Its *path* is γ_0 = −1, γ_j = γ_(j−1) + w_j. The tape
lies in the strip [Lo, Hi − 1] if all γ_j do. The *tape machine* has the moves

    left    … x [s] y …  →  … [x] s y …           (if the head is not on the first letter)
    right   … x [s] y …  →  … x s [y] …           (if the head is not on the last letter)
    split   [s]  →  x′ [−] y′   with x′ + y′ = s + 1:      [+] → + [−] +,     [−] → + [−] −,     [−] → − [−] +
    merge   x [+] y  →  [x + y + 1]:      − [+] − → [−],    − [+] + → [+],    + [+] − → [+]   (head needs both neighbors)

and a move is allowed iff the new tape lies in the strip. There is no merge at + [+] + (x + y + 1 = 3 is not a
step). Split and merge do not change the two ends of the path. The *root* is the tape + [−] +.

**Proposition T.** Let Φ(h, L, U) be the tape whose path is

    −1, l_1, …, l_k, u_j − 1, …, u_1 − 1, 0        (L = (l_1, …, l_k), U = (u_1, …, u_j))

with the head on the step from l_k to u_j − 1 (from the initial −1 if L is empty, to the final 0 if U is empty).
Then, on the states reachable from t:
(a) Φ(s) is a tape in the strip; its head letter is − iff both tops are unmarked, and + iff exactly one is marked.
(b) Φ(s) = Φ(s′) iff s and s′ have the same stacks. Then h′ = h, or the head letter is + and {h, h′} = {l, l + 2}
    with l the lower top (*twin states*).
(c) Φ maps the successors of s bijectively onto the successors of Φ(s) in the tape machine.
Consequently M_m^[Lo,Hi] is the number of walks of length 2m − 2 from the root to the root in the tape machine on
the strip [Lo, Hi − 1]. Without a cap (strip = all integers) this number is the number of all meanders of order m.

*Proof.* Invariants, by induction over the moves. (I1) Neighbors in a stack differ by 1; a non-empty lower stack
has top h or h − 2, a non-empty upper stack has top h or h + 2. Indeed, before a move to h′ the lower top is
h′ − 1 (for d = +1, the unmarked top h) or one of h′ + 1, h′ − 1 (for d = −1, the top h or h − 2); a pushed label is
h′, which differs by 1 from the old top; a popped lower top is h′ − 1 in both cases (for d = −1 the pop needs it),
and the element under it is h′ − 2 or h′ by the first part. So after the move the lower top is h′ or h′ − 2. The
upper side is the mirror image. (I2) The bottom label of each stack is 0, L can be empty only at h = 1, U only at
h = −1: the bottom 0 of L is popped by a move to h′ = 1; then the lower side is marked, so the next move is d = −1
to h′ = 0, where the lower side can only push 0. Mirror image for U.
So with the sentinels −1 (below L) and 0 (after U), the sequence in the proposition is a nearest-neighbor path:
inside L and inside U by (I1), at the two ends by (I2) (steps +1), and across the head because by (I1) and
"not dead" the pair (lower top l, upper top u) is (h, h), (h − 2, h) or (h, h + 2), where an empty L counts as
l = −1 = h − 2 and an empty U as u = 1 = h + 2. Hence (u − 1) − l = −1 in the first case and +1 in the other two.
This is (a); the strip condition is the condition on the labels. (b): the stacks determine the tape and
conversely; if the head letter is −, h = l; if it is +, h ∈ {l, l + 2}.

(c) Case E: head letter −, l = u = h. Write the path near the head as …, r, h, [−], h − 1, r″, ….
 d = +1 (h + 1 ≤ Hi holds because h = l ≤ Hi − 1). Upper side: no pop (it would need u = h + 2); push h + 1.
 Lower side: push h + 1 if h + 1 ≤ Hi − 1, or pop.
   Lower push: the path becomes …, h, h + 1, [−], h, h − 1, …: the split [−] → + [−] −; it needs the level
   h + 1 ≤ Hi − 1, which is the strip condition for the new tape.
   Lower pop: …, r, [·], h, h − 1, …: the old lower top h is now the path level of the new upper label h + 1;
   the head has moved one letter to the left.
 d = −1 (h − 1 ≥ Lo holds because h = u ≥ Lo + 1). Lower side: no pop; push h − 1. Upper side: push h − 1 if
 h − 1 ≥ Lo + 1, or pop.
   Upper push: …, h, h − 1, [−], h − 2, h − 1, …: the split [−] → − [−] +, with the strip condition h − 2 ≥ Lo.
   Upper pop: the head moves one letter to the right.
Case O: head letter +, u = l + 2; put c = l + 1. Both twins move to h′ = c (the twin h = l has the upper top
marked, so d = +1; the twin h = l + 2 has the lower top marked, so d = −1), and c is in the band. On the lower
side push c (allowed: c = u − 1 ≤ Hi − 1) or pop (allowed in both twins: the top is l = h′ − 1) unless L is empty;
on the upper side push c (allowed: c ≥ Lo + 1) or pop (the top is u = h′ + 1) unless U is empty. With the path
…, r, l, [+], l + 1, r″, …:
   push, push: …, l, l + 1, [−], l, l + 1, …: the split [+] → + [−] +.
   push, pop: head one letter to the right.    pop, push: head one letter to the left.
   pop, pop: …, r, [·], r″, …: the three steps (r → l), (l → l + 1), (l + 1 → r″) are replaced by the one step
   r → r″. The new state is dead iff both new tops are marked, iff r = c − 2 and the new upper top is c + 2, iff
   r″ − r = 3, iff the three steps were + + +. Otherwise r″ − r = ±1 and this is the merge.
An empty stack has no pop, which is "the head cannot leave the tape" and "the merge needs both neighbors".
In every case the listed results are all successors, distinct moves give distinct tapes, and the list is exactly
the list of tape moves allowed in the strip. Φ(t) is the path −1, 0, −1, 0, the root.
Count: by (c), a walk from t lifts step by step and uniquely from a tape walk from the root, and it ends at t iff
the tape walk ends at the root (the root has head letter −, so h = 0 is determined). All meanders of order m lie in
the band [−2m, 2m]. ∎

Checks [data, `machine.py`, `machine.log`]: label machine and tape machine, both written from the definitions,
give the same counts for m ≤ 8 for the bands [−1,1], [−1,2], [−2,2], [−3,3] and without cap (Baxter; 1, 2, 7, 31,
…, 38891; 1, 2, 8, 40, …, 78726; …, 105272; …, 110954), and Φ maps successor sets onto successor sets at every
state within 12 moves (40,976 states and 4,970 twin pairs for rung 2; 445,067 states without cap; 0 mismatches).

**Corollary T1 (twins).** Twin states have the same successors. So a function V with Σ_(s′) V(s′) = λ V(s), λ ≠ 0,
takes the same value on twins: eigenfunctions of the label machine are functions on tapes. In rung 2:
V(h = 2; 0*, 2) = V(h = 0; 0, 2*), V(h = 1; −1*, 1) = V(h = −1; −1, 1*), V(h = 0; −2*, 0) = V(h = −2; −2, 0*).

**Dictionary.** Level γ of the path = lower label = upper label − 1. For rung 2 the strip is {−2, −1, 0, 1}; I
number its levels 1, 2, 3, 4 (level = label + 3). The part of the tape left of the head is the lower stack, the
part right of the head is the upper stack read from the top. The dead state is the head on the step 2 → 3 between
a step 1 → 2 and a step 3 → 4: the lower arc under the top has label −2 and the upper arc under the top has label 2.
A strip with fewer than four levels contains no + + +: this is "dead states exist iff Hi − Lo ≥ 4".

## 2. The eigen-equation and the bound lemma [proved]

A *strip tape* is a tape whose path may start at any level of the strip (and stays in it). For a function V on
strip tapes the bulk eigen-equation at a tape T = … x [s] y … whose head has both neighbors is

    (EE)   λ V(T) = V(left T) + V(right T) + Σ over the allowed splits V + V(merge T),

the last term being absent if s = − or x = s = y = +.

**Lemma B.** Let V ≥ 0 on the strip tapes of the strip [Lo, Hi − 1] with Σ_(moves) V ≤ λ V at every tape whose head
has both neighbors (a *supersolution*), and V(− + [−] + −) > 0 for the tape with path 0, −1, 0, −1, 0, −1. Then
M_m^[Lo,Hi] ≤ λ^(2m−2) for all m ≥ 1.

*Proof.* For a tape T of the tape machine (path from −1 to 0) let T⁺ be T with a letter − put in front (the path
now starts at level 0) and a letter − put behind (it ends at −1); T⁺ is a strip tape because Lo ≤ −1 ≤ 0 ≤ Hi − 1.
Its head has both neighbors, and every move of T is a move of T⁺ with result (T′)⁺. So W(T) := V(T⁺) satisfies
Σ_(T → T′) W(T′) ≤ Σ over all moves of T⁺ V ≤ λ W(T). Iterating, (number of walks root → root of length n)·W(root)
≤ λ^n W(root), and W(root) > 0. Proposition T. ∎

## 3. Matrix-product weights [proved]

Let e run over the directed edges (p → q), |p − q| = 1, of the levels of the strip. A *matrix-product weight* is

    V(e_1 … e_k [s] f_1 … f_j) = ⟨W_(p0)| A(e_1) ⋯ A(e_k) · C(s) · B(f_1) ⋯ B(f_j) |W′_(q0)⟩,

with row vectors ⟨W_p| (p = first level of the path), column vectors |W′_q⟩ (q = last level), and matrices A(e),
B(e), C(e). The matrices may be *graded*: a space 𝓛_p for every level on the left and ℛ_q on the right, A(p → q)
of size dim 𝓛_p × dim 𝓛_q, B(p → q) of size dim ℛ_p × dim ℛ_q, C(p → q) of size dim 𝓛_p × dim ℛ_q. Graded matrices
are a special case of D × D matrices with D = Σ dim (blocks). D = 1 is the product weight of 622.

For a walk x s y of three edges write I_s for the split terms without their outer factors,

    s = (a → a − 1):   I_s = [a + 1 in the strip] A(a, a+1) C(a+1, a) B(a, a−1) + [a − 2 in the strip] A(a, a−1) C(a−1, a−2) B(a−2, a−1),
    s = (a → a + 1):   I_s = A(a, a+1) C(a+1, a) B(a, a+1),

and C_(m(x,y)) for the merge: for s = (a → a + 1), x = (r → a), y = (a + 1 → r″): C(r → r″) if |r − r″| = 1 and 0
if r″ − r = 3; for s down, 0. The relations are, one for every walk x s y of three edges in the strip:

    (R)    λ A_x C_s B_y = C_x B_s B_y + A_x A_s C_y + A_x I_s B_y + C_(m(x,y)).

The four terms on the right are the moves left, right, split, merge.

**Proposition S (sufficiency).** If (R) holds, then for all boundary vectors the matrix-product weight satisfies (EE)
at every strip tape whose head has both neighbors.
*Proof.* Put the relation for (x, s, y) between ⟨W| A(e_1) ⋯ A(e_(k−1)) and B(f_2) ⋯ B(f_j) |W′⟩. ∎

**Proposition F (finite rank = matrix product).** For a function V on strip tapes let 𝓛_p be the span of the
functions L ↦ V(L w [s] R) on the left words L ending at level p (w, s, R arbitrary), and ℛ_q the span of the
functions R ↦ V(L [s] w R) on the right words starting at level q. If all these spaces have finite dimension
(*finite rank*), then V is a graded matrix-product weight with these spaces, and V satisfies (EE) at all tapes
whose head has both neighbors if and only if its matrices satisfy (R).
*Proof.* The standard construction for recognizable series. Let α(L) ∈ 𝓛_p* be the evaluation at L. The map
φ ↦ φ(· e) sends 𝓛_q to 𝓛_p for e = (p → q) (by the definition of the spans); let A(e) be its transpose, so
α(L e) = α(L) A(e). Same on the right with evaluations β(R) and matrices B(e). For fixed R the function
V(· [s] R) lies in 𝓛_p, and its coordinates, as functions of R, lie in ℛ_q; so V(L [s] R) = α(L) C(s) β(R) with
C(s) ∈ 𝓛_p ⊗ ℛ_q. With one-edge and empty words the boundary vectors are the evaluations at the empty words.
(EE) at the tape L x [s] y R is α(L)·(λ A_x C_s B_y − right side of (R))·β(R) = 0; the evaluations α(L) span 𝓛_r*
(r the first level of x) and the β(R) span the dual of ℛ at the last level of y, so (EE) at all such tapes is the
matrix identity (R). ∎

So the question "is there a matrix-product eigenfunction with a finite auxiliary space" is exactly the question
"is there a solution of (R) in graded matrices of some size". Two separate chains with a coupling at the tops, one
chain through both stacks, weights on pairs of consecutive labels, auxiliary spaces indexed by the label: all are
instances.

## 4. Width 3 (band [−1, 2]): the D = 1 solutions of (R) are exactly the family of 622 [proved]

Three levels 1, 2, 3 (labels −1, 0, 1 below; 0, 1, 2 above), edges 12, 21, 23, 32, eight walks x s y. D = 1, all
weights nonzero. Gauge-invariant letter weights: p_− = a_21 a_12 and p_+ = a_23 a_32 (a visit of the lower stack to
label −1 or 1: the letters − and + of σ in 622), q_+ = b_21 b_12 and q_− = b_23 b_32 (a visit of the upper stack to
label 0 or 2: the letters + and − of τ).
The relations with s = 21 (x = 12 or 32; y = 12) are λ a_x c_21 b_12 = c_x b_21 b_12 + a_x a_21 c_12 + a_x a_23 c_32 b_21 b_12;
their difference gives c_12/a_12 = c_32/a_32 =: P. The relations with s = 32 (x = 23; y = 21 or 23) give in the same
way c_21/b_21 = c_23/b_23 =: Q. With these, the two relations of each s coincide (for s = 12 the merge gives c_y,
for s = 23 it gives c_x) and, after division by b_21, a_32, a_12, b_23, the four remaining relations are

    s = 21:  λ Q = P (1 + p_+ + p_−/q_+)          s = 23:  λ Q = P (1 + p_+)(1 + 1/q_−)
    s = 32:  λ P = Q (1 + q_+ + q_−/p_+)          s = 12:  λ P = Q (1 + q_+)(1 + 1/p_−).

The first line gives p_− q_− = q_+ (1 + p_+), the second p_− q_− = p_+ (1 + q_+). Hence

    p_+ = q_+ =: p,     p_− q_− = p (1 + p),     λ² = (1 + p)(1 + p_−)(1 + q_−)/p,     (Q/P)² = p_− (1 + q_−)/(q_− (1 + p_−)),

which is Proposition 3 of 622 with z = 1/λ and κ = Q/P. So the method reproduces the refereed family, exactly, as
its D = 1 solution set, and nothing else at D = 1.
Check [exact, `width3_exact.py`]: a Gröbner basis of the eight relations contains (p_− + 1)(p_+ − q_+)(q_+ + 1),
(p_− + 1)(p_− q_− − q_+² − q_+) and λ² p_+ − (p_− q_− + p_− q_+ + p_− + p_+² q_+ + p_+² + p_+ q_− + q_− + q_+ + 1); on the family all
eight relations vanish identically.

## 5. Rung 2 at D = 1: the relaxed family and Λ_2 ≤ 6 + 4√2 [proved]

Four levels, edges 12, 21, 23, 32, 34, 43, sixteen walks x s y. Call the walk (12, 23, 34) the *dead context*.
As in section 4, the relations with s down force c_12/a_12 = c_32/a_32 =: P_2 and c_23/a_23 = c_43/a_43 =: P_3 (from
s = 21 and s = 32), c_21/b_21 = c_23/b_23 =: Q_2 and c_32/b_32 = c_34/b_34 =: Q_3 (from s = 32 and s = 43). Then for
s = 23 the merge terms of the three live contexts are equal, c_m/(a_x b_y) = P_2/b_32, and the dead context has 0
in its place. So the sixteen relations cannot all hold (this is Proposition 4 of 622), but fifteen can, and then
the sixteenth holds with "<": a supersolution.

**Proposition 5.** The positive D = 1 weight (gauge a_12 = a_23 = a_34 = 1, b_21 = b_32 = b_43 = 1)

    a_21 = b_34 = 1,   a_32 = b_23 = 1/√2,   a_43 = b_12 = √2 − 1,
    c_12 = 1,  c_32 = c_34 = 1/√2,  c_23 = ρ,  c_43 = ρ(√2 − 1),  c_21 = ρ√2,     ρ = (2 + √2)/4,     λ = 2 + √2

satisfies the fifteen live relations with equality and the dead context with ratio (right side)/(left side) =
4√2 − 5 = 0.6569. Hence, by Lemma B,

    M_m^(2) ≤ (6 + 4√2)^(m−1) = 11.6568542…^(m−1)    for all m ≥ 1.

*Proof.* With A_1 = a_21, A_2 = a_32, A_3 = a_43, B_1 = b_12, B_2 = b_23, B_3 = b_34, P_2 = 1, P_3 = ρ (so Q_3 = A_2,
Q_2 = ρ/B_2) the six equations of the live contexts are

    (21)  λ ρ/B_2 = 1 + A_1/B_1 + A_2                    (12)  λ = (ρ/B_2)(1 + B_1)(1 + 1/A_1)
    (32)  λ A_2 = ρ(1 + A_3) + (ρ A_2/B_2)(1 + B_1)      (23)  λ ρ = (1 + A_2)(1 + B_2)
    (43)  λ ρ = A_2 B_3/A_3 + A_2 (1 + B_2)              (34)  λ A_2 B_3 = ρ (1 + A_3)(1 + B_3).

Put A_1 = B_3 = u, A_2 = B_2 = v, A_3 = B_1 = w. (23): λρ = (1 + v)². (43) and (21) both become u/w = (1 + v)/v.
(12) and (34) both become λ = ρ(1 + w)(1 + u)/(u v). (32) is λ = 2ρ(1 + w)/v. So u = 1, w = v/(1 + v), and
λ² = 2(1 + w)(1 + v)²/v = 2(1 + 2v)(1 + v)/v = 2(1/v + 3 + 2v), smallest at v = 1/√2 with value 6 + 4√2. All
weights are positive. The dead context: by the first paragraph its right side is the common value minus the
merge term. Lemma B applies (the weight is a product of positive numbers, so V > 0). ∎

Checks [exact, `d1_bound.py`, `d1_bound.log`]: the sixteen ratios in Q(√2) (fifteen equal 1, one equals 4√2 − 5);
sub-invariance at all 17,812 tapes within 12 moves of the root (largest ratio 1.000000000000001, floating point);
the 27 stored terms satisfy M_m ≤ (6 + 4√2)^(m−1) in integer arithmetic (M_27/(6 + 4√2)^26 = 2.7·10⁻⁵).
Remarks. 606 and 622 quote this radius-1 value as "(1 + √2)² = 11.657"; the number 11.657 is (2 + √2)² = 6 + 4√2 =
2(1 + √2)². No cone angle is attached to it: the weight is not an eigenfunction, so the ratio of two members of
the relaxed family is not a martingale of any chain of the machine.
The same computation without levels (all meanders, weights a_±, b_±, c_±) gives λ = t + 1/t + 2 ≥ 4, the trivial
bound 16.

## 6. The no-go theorem for matrix products with cancellable push matrices [proved]

**Theorem N.** Let the strip contain four consecutive levels, called 1, 2, 3, 4 (so the band has width at least 4;
all meanders are included). Let graded or ungraded matrices satisfy (R) with some λ, and assume
 (H1) the square matrix B_21 B_12 is invertible, and
 (H2) A_12 is left-cancellable (A_12 X = 0 implies X = 0; for a square matrix: invertible).
Then C_34 = 0: the weight vanishes at every tape with the head on a step 3 → 4.
Mirror statement: if A_34 A_43 is invertible and B_34 is right-cancellable, then C_12 = 0.
In particular no solution of (R) with invertible matrices A(e), B(e), of any size D, gives a positive weight.
By Proposition F this covers every matrix-product weight with invertible matrices whose left vectors and right
vectors span the auxiliary spaces.

*Proof.* The walks (x, 21, 12) with x = 12 and x = 32 give, by (R) (no merge, s is down),

    C_x B_21 B_12 = A_x Z,     Z := (λ C_21 − I_21) B_12 − A_21 C_12   (the same Z for both x).

By (H1), C_x = A_x P for x = 12 and x = 32, with P = Z (B_21 B_12)^(−1).
The walks (x, 23, y) with x ∈ {12, 32} and y ∈ {32, 34} give, by (R) and C_x = A_x P,

    C_(m(x,y)) = A_x Γ_y,     Γ_y := λ C_23 B_y − P B_23 B_y − A_23 C_y − I_23 B_y,

with C_(m(32,32)) = C_32, C_(m(32,34)) = C_34, C_(m(12,32)) = C_12 and C_(m(12,34)) = 0 (the dead context).
The dead context says A_12 Γ_34 = 0, so Γ_34 = 0 by (H2), so C_34 = A_32 Γ_34 = 0. ∎

Remarks.
1. For D = 1 and nonzero weights (H1), (H2) hold, and the theorem is Proposition 4 of 622 (positivity is not
   needed there, only nonvanishing).
2. A finite linear combination of product weights with nonzero letter weights is the case of diagonal
   invertible matrices A(e), B(e) (after merging proportional terms; C may then be rectangular). So no such
   combination (a "Bethe sum" over letter weights) is a positive eigenfunction.
3. By Proposition F, a positive eigenfunction of finite rank, in its minimal realization, must violate (H1) or
   (H2) and their mirror images: either the left vectors α(L·(1 → 2)) of the words that end by coming up from the
   bottom level do not span 𝓛_2* (so the left environment cannot be a product; it has to carry the information
   "the label under the top is the lowest one" in a way that restricts the whole vector), or the double step
   2 → 1 → 2 on the right is a singular operator. Weights with finite memory (factors that depend on the last k
   steps) are of this kind; they are not excluded by the theorem. Section 4 of FINAL.md has the numerical search
   in that class.

## 6b. Weights with memory one collapse to product weights [proved]

The smallest singular realizations are the weights with *memory one*: one factor per pair of consecutive label
edges on each side, and a head weight that sees its two neighbors,

    V(e_1 … e_k [s] f_1 … f_j) = ω(e_1) · Π_(i<k) f(e_i, e_(i+1)) · c(e_k, s, f_1) · Π_(i<j) g(f_i, f_(i+1)) · ω′(f_j),     k, j ≥ 1,

all numbers positive. (In matrix form: the left space at a level has one basis vector per edge into the level.
This contains the "product times window of radius 2" of 622 and the second-order label factors of 606.)

**Theorem N2.** On the strip with four levels (rung 2) no positive weight with memory one satisfies (EE) at all
tapes with k ≥ 2 and j ≥ 2. More precisely, (EE) at these tapes forces f(X, e) = φ(X)ψ(e), g(e, Y) = γ(e)δ(Y) and
c(x, s, y) = c_s φ(x) δ(y): the weight is a product weight, and Theorem N (D = 1) applies.

*Proof.* Divide (EE) at the tape … x′ x [s] y y′ … by the common factors. With U(x′, x, s) := c(x′, x, s)/f(x′, x),
W(s, y, y′) := c(s, y, y′)/g(y, y′) and S(x, s, y) := λ c(x, s, y) − Σ_(splits) f(x, i_1) c(i_1, −, i_2) g(i_2, y):

    (★)   S(x, s, y) = U(x′, x, s) g(s, y) + f(x, s) W(s, y, y′) + c(x′, m, y′)/(f(x′, x) g(y, y′)),

the last term only if the merge exists; m is the merged step, so the last term is the head weight of the tape
… x′ [m] y′ …. The left side does not depend on (x′, y′).

Step 1 (s down, no merge). U(x′, x, s) does not depend on x′, and W(s, y, y′) does not depend on y′. In words:
(D-L) c(X, H, Y)/f(X, H) does not depend on X if Y is down; (D-R) c(X, H, Y)/g(H, Y) does not depend on Y if X is down.

Step 2 (head 32). Put M(X, Y) := c(X, 32, Y)/(f(X, 32) g(32, Y)) for X ∈ {23, 43}, Y ∈ {21, 23}. (D-L) with Y = 21 gives
M(23, 21) = M(43, 21); (D-R) with X = 43 gives M(43, 21) = M(43, 23). (★) at … x′ (32) [23] (32) y′ …, where the merge
gives the head 32 between x′ and y′, reads S(32, 23, 32) = U(x′, 32, 23) g(23, 32) + f(32, 23) W(23, 32, y′) + M(x′, y′): a
function of x′ plus a function of y′ plus M is constant, so the second difference of M vanishes, so
M(23, 23) = M(23, 21). Hence
    (P32)   c(X, 32, Y) = M_0 f(X, 32) g(32, Y)   for all four (X, Y),
and U(x′, 32, 23) = M_0 g(32, 23), W(23, 32, y′) = M_0 f(23, 32) do not depend on x′, y′.

Step 3. (★) at … (21)(12) [23] (32) y′ …, y′ ∈ {21, 23} (merge + [+] − → head 12 between 21 and y′): the first two terms
are constant by Step 2, so c(21, 12, y′)/g(32, y′) does not depend on y′. By (D-R), c(21, 12, y′)/g(12, y′) does not
depend on y′. Hence (G2): g(12, Y)/g(32, Y) does not depend on Y ∈ {21, 23}.
Mirror image (reverse the tape and exchange the levels p ↔ 5 − p; this exchanges left and right, 12 ↔ 34, 21 ↔ 43,
fixes 23 and 32, and maps the tape machine to itself): (F3): f(X, 34)/f(X, 32) does not depend on X ∈ {23, 43}.

Step 4 (head 23). (★) at … x′ (21) [12] (21)(12) …, x′ ∈ {12, 32} (merge − [+] − → head 21 between x′ and 12):
S(21, 12, 21) = U(x′, 21, 12)·(g(12, 21) + 1/g(21, 12)) + f(21, 12) W(12, 21, 12), so U(x′, 21, 12) does not depend on x′.
(★) at … x′ (21) [12] (23) y′ …, x′ ∈ {12, 32}, y′ ∈ {32, 34} (merge − [+] + → head 23 between x′ and y′; all four tapes
are live):
    S(21, 12, 23) = U(x′, 21, 12) g(12, 23) + f(21, 12) W(12, 23, y′) + N(x′, y′),     N(x′, y′) := c(x′, 23, y′)/(f(x′, 21) g(23, y′)).
The first term is constant, so N(x′, y′) =: n(y′) does not depend on x′; and W(12, 23, y′) = c(12, 23, y′)/g(23, y′) =
n(y′) f(12, 21), so (1 + f(21, 12) f(12, 21))·n(y′) is constant: n(32) = n(34) =: n_0. Hence
    (P23a)   c(X, 23, Y) = n_0 f(X, 21) g(23, Y)    for X ∈ {12, 32}, Y ∈ {32, 34},
and by the mirror image
    (P23b)   c(X, 23, Y) = n_0′ f(X, 23) g(43, Y).
Comparing, f(X, 21)/f(X, 23) = (n_0′/n_0)·g(43, Y)/g(23, Y); the left side depends only on X, the right side only on
Y, so both are constant: (F2) f(X, 21)/f(X, 23) does not depend on X ∈ {12, 32}; (G3) g(43, Y)/g(23, Y) does not
depend on Y ∈ {32, 34}.

Step 5. Levels 1 and 4 have one edge in and one edge out. By (F2), (F3) the pair factors on the left are
f(X, e) = φ(X)ψ(e) at every level, and by (G2), (G3) g(e, Y) = γ(e)δ(Y). So the left part of the weight is a
boundary factor of e_1 times Π a(e_i)/φ(e_k) with a = φψ, the right part is Π b(f_i)/δ(f_1) times a boundary factor
of f_j, with b = γδ, and the head carries c̃(x, s, y) = c(x, s, y)/(φ(x)δ(y)). This c̃ does not depend on (x, y):
for s = 32 by (P32); for s = 23 by (P23a); for s = 21 because U(x′, 21, 12) is constant (Step 4; y = 12 is the only
choice); for s = 12 because c(21, 12, Y)/g(12, Y) is constant ((D-R); x = 21 is the only choice); s = 43 and s = 34
are the mirror images. So on the tapes with k, j ≥ 1 the weight is a D = 1 product weight times a factor that
depends only on the first and the last letter of the tape, and that factor is common to a tape with k, j ≥ 2
and all the tapes its moves lead to. The D = 1 relations (R) follow, and Theorem N at D = 1 gives c_34 = 0, a
contradiction with positivity. ∎

Check [data, `memory1_check.py`]: least squares on the 42 equations (★) with positive unknowns has floor 2.6·10⁻² at
λ = 3.38 (6.2·10⁻³ at λ = 4.0) on four levels; on three levels the same program finds exact solutions (8·10⁻¹⁷) for
λ ≥ 3.2288, with pair factors of rank one, and the D = 1 floor 8.9·10⁻³ at λ = 3.20.

Remarks. The proof is a collapse: memory one implies memory zero. The same scheme suggests the statement
"memory k implies memory k − 1" for every k, which would exclude every weight with finite memory; I did not
attempt k = 2. The four-level structure is used (levels 1 and 4 have a single neighbor).

## 7. Supersolutions in matrix-product form: the certificate statement [proved]

**Proposition C.** Let A(e), B(e), C(e) be D × D matrices with positive entries, indexed by the six edges of the four
levels of rung 2, and let

    λ* = max over the sixteen walks x s y and all entries (i, j) of  (right side of (R) without λ)_(ij) / (A_x C_s B_y)_(ij).

Then M_m^(2) ≤ (λ*)^(2m−2) for all m.
*Proof.* With boundary vectors with positive entries the weight is positive, and Σ_(moves) V − λ* V at a tape
L x [s] y R is α(L)·(right side − λ* A_x C_s B_y)·β(R) with α(L), β(R) ≥ 0 entrywise and a matrix ≤ 0 entrywise.
Lemma B. ∎
