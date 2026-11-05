# 630 REPAIRS: the arguments that the chain left thin, written out in full (2026-10-10)

Referee notes for the chain 204 (Lemmas 1 to 3) → 606 section 1 → 622 sections 0, 2, 3 → 625 sections 0 to 4.
No step of the chain was found to be wrong. Each section below supplies an argument that the source gives only in
compressed form, only by a script, or not at all. Status of every item here: **proved** (complete argument in this
file), unless marked otherwise. Scripts named here are in this folder and were written from the definitions; they
are checks, not parts of the proofs.

## 0. Conventions

- A *system* of order m is a pair (u, l) of noncrossing perfect matchings of the points 1, …, 2m on a line; u is drawn
  above the line (upper arcs), l below (lower arcs). A *meander* is a system whose arcs form one closed curve.
  M_m^all = 1, 2, 8, 42, 262, … is the number of meanders of order m.
- Point 2m is called ∞. An arc is *finite* if it does not end at 2m. The upper arc at 2m is (p, 2m), the lower one
  (q, 2m); they are called the *rays*. The river edges e_i = (i, i+1), 1 ≤ i ≤ 2m − 2, are the *finite edges*.
- For a meander, the finite arcs form a path through all the points 1, …, 2m − 1 (the closed curve with the point 2m
  removed). The *winding function* g is the unique function on {1, …, 2m − 1} with g(1) = 0, g(b) = g(a) − 1 for every
  finite upper arc a < b and g(b) = g(a) + 1 for every finite lower arc a < b. It exists and is unique because the
  finite arcs form a path. No property of g is assumed in this definition.
- The band family: meanders with −1 ≤ g(i) ≤ 2 for all i ≤ 2m − 1. M_m is their number. It depends on the choice of
  the point 2m as the root: g is not invariant under rotation of the root (204 section 6). The family is a set of
  rooted objects in the same sense in which M_m^all counts rooted meanders.
- U means a step +1 of g along the river, D a step −1.
- The *face of a finite upper arc* (a, b) is the region between the arc, the river and the upper arcs directly nested
  in (a, b) (its *children*). Its *floor* is the set of river edges on its boundary: if the children are
  (a_1, b_1), …, (a_k, b_k) from left to right, the floor edges are (a, a_1), (b_1, a_2), …, (b_k, b); there are k + 1
  of them (for k = 0 the single edge (a, b)). The two *unbounded upper faces* lie left and right of the upper ray; the
  finite floor edges of the left one are the edges (b_i, a_(i+1)) between consecutive top-level arcs left of p and the
  edge from the last of them to p; if there are j such arcs there are j such edges. The same on the right (j′ arcs,
  j′ edges). The same definitions below the river.

## 1. Link 1: the winding bridge (204, Lemmas 1 to 3)

Verdict: correct. Two places of the proof of Lemma 3 are compressed; they are written out here. For completeness the
statement that the later links use is isolated as Theorem W.

**Theorem W.** (a) For a meander, g has steps ±1 along the river, g(2m − 1) = 0 and g(p) = g(q) = 0; every bounded
upper face has exactly one D among its floor edges and the two unbounded upper faces have none among their finite
floor edges; every bounded lower face has exactly one U and the unbounded lower faces none.
(b) Conversely, if (u, l) is a system and G is a function on {1, …, 2m − 1} with G(1) = 0, steps ±1 along the river,
and the two arc rules, then (u, l) is a meander and G = g.

*Proof of (a).* This is 204 Lemma 1 followed by Lemma 2; I checked both line by line and have nothing to add except
the following remarks, which make the argument self-contained.
- The road is a simple closed curve through ∞ on the sphere; it cuts the sphere into two discs. The river crosses
  the road at each of the 2m points, so each river edge (including the two that end at ∞) lies in one disc and is a
  chord of it; chords in one disc are disjoint, so the faces in a disc and its chords form a tree. The two river
  edges at ∞ lie in different discs because the river crosses the road at ∞. Removing the chord at ∞ from the tree
  of a disc leaves two trees, and the two faces of that disc next to the removed chord are exactly the two faces of
  the disc that have ∞ on their boundary. These four faces are the two unbounded upper and the two unbounded lower
  faces. Root the four trees there. Every finite edge e then has a child face (the face on the side away from the
  root), and it is either the face above e or the face below e. Put s_e = D in the first case and U in the second.
- A bounded upper face F is a non-root vertex of its tree, so it is the child of exactly one finite edge, and that
  edge is on its floor; an edge on the floor of F is D if and only if F is its child. So F has exactly one D on its
  floor. An unbounded upper face is a root and is the child of no edge. The same below with U. (Lemma 1.)
- Let G be the walk with G(1) = 0 and steps s. For a finite upper arc (a, b) with children (a_i, b_i),
  G(b) − G(a) = (sum of the steps on the floor) + Σ_i (G(b_i) − G(a_i)). By induction from the innermost arcs,
  G(b_i) − G(a_i) = −1, and the floor has k + 1 steps of which exactly one is D, so G(b) − G(a) = (k − 1) − k = −1.
  The same below gives +1. So G obeys both arc rules and G(1) = 0, hence G = g by uniqueness. The top-level steps
  are all U (they are the finite floor edges of the unbounded upper faces), so
  g(p) − g(1) = j − j = 0 and g(2m − 1) − g(p) = j′ − j′ = 0. Below the river the top-level steps are all D and the
  lower arcs raise g by one, which gives g(q) = 0 in the same way. ∎

*Proof of (b).* Step 1 (one marked edge per face). For a finite upper arc with k children, the arc rule for it and
for its children gives (sum of the k + 1 floor steps) = k − 1, so exactly one floor step is D. So each of the m − 1
bounded upper faces has exactly one D on its floor, and likewise each of the m − 1 bounded lower faces exactly one
U. Floors of different upper faces are disjoint, and so are floors of different lower faces. Hence there are at
least m − 1 steps D and at least m − 1 steps U among the 2m − 2 finite edges; so there are exactly m − 1 of each,
every D is the unique D of a bounded upper face, every U is the unique U of a bounded lower face, and the finite
floor edges of the unbounded upper faces carry no D and those of the unbounded lower faces no U. Call the face
selected by a finite edge (the upper face above it if the edge is D, the lower face below it if it is U) its
*target*. Every bounded face is the target of exactly one edge, and that edge is on its boundary.

Step 2 (no second loop). Suppose the arcs form k ≥ 2 closed curves. They are disjoint simple closed curves on the
sphere, so the complement has k + 1 regions, and the graph with the regions as vertices and the curves as edges is a
tree with k edges. Let C_∞ be the curve through ∞. A tree has at least two leaves. If every leaf had C_∞ as its
only boundary curve, the two ends of the edge C_∞ would both be leaves and the tree would be a single edge, so
k = 1. Hence there is a leaf region R whose boundary is one curve C ≠ C_∞. R is an open disc and ∞ is not in its
closure. A river edge has no point of a curve in its interior, so it lies in one region; the river edges inside R
are disjoint chords of R with both ends on C, they are finite edges (an edge ending at ∞ has an end on C_∞), and
there are no other arcs or river edges inside R. If there are r such chords, they cut R into r + 1 faces of the
system. Each of these faces is bounded (it does not have ∞ on its boundary), so by Step 1 each is the target of
exactly one finite edge on its boundary; the river edges on the boundary of a face inside R are among the r
chords, and an edge has one target. So r + 1 faces need r + 1 different chords. Contradiction. So k = 1.

Step 3. (u, l) is a meander and G obeys the defining rules of g, so G = g. ∎

Checks (`link1_check.py`, `brute.c`; logs with the same stems). For every system of order m ≤ 7 and every ±1 walk
(not assumed to be a bridge): the systems that admit a walk obeying both arc rules are exactly the meanders, each
with exactly one walk, equal to g computed along the road (13,820 meanders at m = 7). For every meander of order
m ≤ 11 (73,424,650 at m = 11): g computed along the road has ±1 river steps, g(2m − 1) = 0, g(p) = g(q) = 0; no
violation. For m ≤ 7 the rule of Lemma 1 in the form "the step at e_i is D if and only if the part of the road
between the points i and i + 1 that avoids ∞ leaves i by its upper arc" reproduces the steps of g.

What this link does not cover: Lemma 4 of 204 and sections 2 to 8 there were not refereed.

## 2. Link 2: the label machine (606 section 1, in the label form of 622 section 0)

Verdict: correct after repair. The gap: ledger row AE1 and the Claim of 606 are stated for the symmetric rungs
max |g| ≤ W, in the form with one bit per open arc; 622 section 0 restates the machine with explicit labels for an
arbitrary band [Lo, Hi] and cites AE1 for it. The band form was checked there by counts only. The argument (a) to
(d) of 606 is correct in substance and carries over, but it was never written for bands or for the label form. Here
is a direct proof.

**The machine.** Fix integers Lo ≤ −1 and Hi ≥ 1. A state is (h, L, U): h an integer in [Lo, Hi], L and U finite
lists of integers (the *lower stack* and the *upper stack*, bottom first). t = (0, (0), (0)). A *move* from (h, L, U)
to (h′, L′, U′):
- h′ = h + d with d = ±1 and Lo ≤ h′ ≤ Hi;
- (step condition) if d = +1, L is not empty and its top is h; if d = −1, U is not empty and its top is h;
- L′ is L with h′ pushed, allowed if h′ ≤ Hi − 1, or L with its top popped, allowed if L is not empty and its top
  is h′ − 1;
- U′ is U with h′ pushed, allowed if h′ ≥ Lo + 1, or U with its top popped, allowed if U is not empty and its top
  is h′ + 1.

(622 also forbids entering a state with no successor. This changes no walk that ends at t.)

**Theorem M.** For a meander (u, l) of order m with Lo ≤ g ≤ Hi put, for 1 ≤ j ≤ 2m − 1,

    S_j = (g(j), L_j, U_j),   L_j = the list of g(a) over the lower arcs (a, b) with a ≤ j < b, in increasing order of a,

and U_j the same for the upper arcs (the rays count as arcs with b = 2m). Then Φ(u, l) = (S_1, …, S_(2m−1)) is a walk
of the machine from t to t with 2m − 2 moves, and Φ is a bijection from the meanders of order m in the band onto
the set of such walks. In particular M^[Lo,Hi]_m is the number of walks of 2m − 2 moves from t to t.

*Proof.* (1) *Φ(u, l) is a walk from t to t.* The arcs open after a point are nested, and the order by left end is
from the outermost to the innermost, so the last entry of L_j belongs to the innermost open lower arc.
S_1 = t: the arcs at point 1 are open after point 1 and g(1) = 0. S_(2m−1) = t: after the point 2m − 1 only the rays
are open, and g(2m − 1) = g(p) = g(q) = 0 by Theorem W (a).
Let 1 ≤ j ≤ 2m − 2, h = g(j), h′ = g(j + 1), d = h′ − h; d = ±1 by Theorem W (a), and h′ is in the band.
*Stacks.* On the lower side the point j + 1 is either the left end of a lower arc, and then L_(j+1) is L_j with
g(j + 1) = h′ pushed; or the right end of a lower arc (a, j + 1), which is the innermost arc open after j, and then
L_(j+1) is L_j with its top g(a) popped, and g(a) = h′ − 1 by the arc rule, so the pop is allowed. In the first case
the push is allowed: if the new arc (j + 1, b) is finite then h′ + 1 = g(b) ≤ Hi; if it is the ray then h′ = g(q) = 0
≤ Hi − 1. The same above: a pushed label h′ satisfies h′ − 1 = g(b) ≥ Lo, or h′ = g(p) = 0 ≥ Lo + 1.
*Step condition.* Let d = +1. The edge e_j is on the floor of exactly one lower face: the face of the innermost
lower arc A = (a, b) open after j if there is one, the unbounded lower face left of the lower ray if there is none.
By Theorem W (a) the U at e_j is the unique U of a bounded lower face. So A exists and is finite, and no floor
edge of A left of e_j is U. Let c be the number of children of A that end at or before j; since A is the innermost
arc open after j, the point j is a (if c = 0) or the right end of the c-th child. The floor edges of A left of the
point j are (a, a_1), (b_1, a_2), …, (b_(c−1), a_c): c edges, all D. So g(j) − g(a) = −c + c = 0: the top of L_j is h.
The case d = −1 is the mirror image.
So S_j → S_(j+1) is a move.

(2) *Φ is injective.* The walk gives g (the first components), and it gives u and l: the point j + 1 is a left end on
a side if and only if the stack of that side grows in the move S_j → S_(j+1); a noncrossing matching with a marked
last point is determined by which of the points 2, …, 2m − 1 are left ends.

(3) *Φ is surjective.* Let (S_1, …, S_(2m−1)) be a walk from t to t, S_j = (h_j, L_j, U_j). Put G(j) = h_j: a ±1 walk
with G(1) = 0 and Lo ≤ G ≤ Hi. Build l: the entry of t belongs to point 1; a push at the move into S_(j+1) makes
j + 1 the left end of a new arc; a pop makes j + 1 the right end of the arc whose left end made the popped entry.
At the end one entry is left (L_(2m−1) has length 1); join its point to 2m. A stack discipline gives a noncrossing
perfect matching. Build u in the same way. By construction the entry made at a point a has the value G(a) (the
push writes h′; the entry of t is 0 = G(1)), so L_j is the list of G(a) over the lower arcs open after j, and a pop at
the point b removes G(a) = G(b) − 1 for the arc (a, b): the lower arc rule. Above, a pop removes G(a) = G(b) + 1: the
upper arc rule. By Theorem W (b), (u, l) is a meander and G = g. It lies in the band, and Φ(u, l) is the given walk. ∎

Remarks.
- The proof of (3) does not use the step condition. So the step condition is implied for walks that return to t
  (by (1)); it removes only walks that cannot return. `machine_check.py` confirms this: the number of closed walks
  is the same with and without it (m ≤ 6, five bands).
- The machine does not know which arc is a ray until the end. The end state t says that the two arcs left open have
  label 0 and have not received a step of their face; an empty stack plays the role of the unbounded face on the
  other side of the ray (no step of that kind is possible there).
- In every state of a walk from t to t: the top of L is h or h − 2, the top of U is h or h + 2, and consecutive
  labels in a stack differ by 1 (from (1): g(j) − g(a) is twice the number of U on the floor of A left of j, and that
  number is 0 or 1; a child opened at j + 1 has label g(j + 1) = g(a) ± 1).
- The hypothesis Lo ≤ −1 ≤ 1 ≤ Hi is used only for the push of a ray. The band [−1, 2] satisfies it.
- The bit form of 606 is the same machine: the bit of the top lower arc with label x is (h − x)/2, and the bit of a
  buried arc is (label of its child − x + 1)/2. I checked the move rules of 606 against the rules above case by
  case; they agree.

Checks (`machine_check.py`). For the bands [−1, 2] (m ≤ 8), [−2, 1], [−1, 1], [−2, 2], [−3, 3] (m ≤ 6): for every
meander in the band (brute force), the sequence S_j read from the definition is a legal walk from t to t; every
walk from t to t found by exhaustive search is rebuilt into a system that is a meander in the band with g equal to
the first components; the two sets of walks are equal (38,891 objects for [−1, 2] at m = 8).

## 3. Link 3: the two-stack form and the upper bound (622 sections 0, 2, 3)

Verdict: correct. I re-derived Proposition 2 by hand from the machine rules of section 2 above with Lo = −1, Hi = 2
(all eight cases: level 1 with d = ±1, the return from h = 2; level 0 with d = ±1, the return from h = −1; the two
boundary cases σ empty and τ = (+)), Proposition 3 and the solution of (E1) to (E4), and Theorem A including the
minimization. Three small points:

**3.1 "Exactly" in Proposition 2.** The induction in 622 proves that every reachable state with h ∈ {0, 1} has the
stated form and that the moves are those of the table. That every state of the stated form is reachable is not
proved there; it follows from step (i) of the proof of Theorem 2 in 625 (the building procedure). Neither bound
needs it beyond that step.

**3.2 The boundary, exactly.** The macro states are (level, σ, τ) with τ not empty and bottom letter + (I write
words bottom first here, as 622 does). The only deviations of the true table from the bulk table are:
- level 1, σ empty: only the two moves "go, σ unchanged, τ keeps or receives +" exist;
- level 0, τ = (+): the two stays in which τ loses its top do not exist (in the raw machine the upper stack is empty
  in the middle of this stay, at h = −1, and the push of 0 that follows recreates the bottom arc).
At these states the sum Σ z^(points)·V(s′)/V(s) over the existing moves is, at the critical member,
z(1 + p) = m₀/(1 + m₀) = 0.4785 and z²(1 + m₀) + z(1 + p) = 0.6624. At all other states it is 1 by (E1) to (E4).
No raw state reached from t has no successor (band width 3). `macro_check.py` compares the table with the raw
machine at every macro state reached within 20 points (94,980 states): equal.

**3.3 The upper bound with moves of two lengths.** Let Q be the one-point transition relation of the raw machine
restricted to the states reachable from t, and extend V to the middle states of stays (h = 2 or h = −1) by
V(s) = z·Σ_(s→s′) V(s′). Then z·Σ_(s→s′) V(s′) ≤ V(s) at every reachable raw state: equality at middle states by
definition; at a macro state the left side is z·Σ_gos V + z·V(middle) = Σ_gos z·V + Σ_stays z²·V ≤ V(s) by 3.2.
Since V > 0, the number N_n of raw walks t → t with n moves satisfies N_n·z^n·V(t) ≤ V(t). With n = 2m − 2 and
Theorem M: M_m ≤ z^(−(2m−2)) = Λ^(m−1). Also, every macro path from t to t with n points has weight
Π z^(points)·V(s′)/V(s) = z^n, so M_m = Λ^(m−1)·P_c(t → t in 2m − 2 points) for the sub-Markov kernel
P_c(s, s′) = z^(points)·V(s′)/V(s) (identity (4.1) of 622).

Exact values (`exact_check.py`, arithmetic in Q(m₀)): 4m₀⁴ = 2m₀ + 1; p = 1/(2m₀); z = 1/(2m₀²(1 + m₀));
(E1) = (E2) = 1; 1/z² = (1 + 2m₀)(1 + m₀)² = Λ; 16Λ⁴ − 112Λ³ − 520Λ² − 536Λ + 25 = 0. The quartic is irreducible
over Q and has two real roots, 0.0446860… and Λ = 10.4245356612221488… (`algebra_check.log`).

## 4. Link 4: the lower bound (625 sections 0 to 4)

Verdict: correct. Two statements are proved in 625 by a script only (the finite check of Theorem 1; Proposition 2
(a)). Both have short hand proofs, given here. The proof of Theorem 2 is complete; 4.4 writes out the bookkeeping.

In this section words are written top first, as in 625. β = m₀/(1 + 2m₀). Recall p = 1/(2m₀), z = p/(m₀(1 + m₀)),
1/(1 + p) = 2m₀/(1 + 2m₀), p/(1 + p) = 1/(1 + 2m₀).

**4.1 Theorem 1 of 625 by hand.** Write F for a stack whose letters are independent fair coins and cF for a stack
with top letter c and independent fair coins below. Under π: at level 1, σ is +F with probability 1 − β and −F with
probability β, and τ is F; at level 0 the same with σ and τ exchanged. The move table at the critical member is
symmetric under the exchange of the stacks together with the exchange of the levels, so it is enough to compute
the law after one move on the event "level 1 at time 1". Given the level and the active top at time 0, the move is
independent of all other letters, so after the move each stack is again of the form F or cF, with all letters
independent. The cases, with their probabilities given "level 1 at time 0" or "level 0 at time 0":

| start | move | probability | σ at time 1 | τ at time 1 |
|---|---|---|---|---|
| level 1, top + | stay, σ keeps, τ keeps | (1 − β)·z² | +F | F |
| level 1, top + | stay, σ keeps, τ receives − | (1 − β)·z²m₀ | +F | −F |
| level 1, top + | stay, σ loses, τ keeps | (1 − β)·z²/p | F | F |
| level 1, top + | stay, σ loses, τ receives − | (1 − β)·z²m₀/p | F | −F |
| level 0, top + | go, τ unchanged, σ keeps | (1 − β)·z | F | +F |
| level 0, top + | go, τ unchanged, σ receives + | (1 − β)·zp | +F | +F |
| level 0, top − | go, τ keeps, σ keeps | β·z | F | −F |
| level 0, top − | go, τ keeps, σ receives + | β·zp | +F | −F |
| level 0, top − | go, τ loses, σ keeps | β·z/m₀ | F | F |
| level 0, top − | go, τ loses, σ receives + | β·zp/m₀ | +F | F |

The probabilities add up to (1 − β)/(1 + m₀) + (1 − β)m₀/(1 + m₀) + β = 1, so the level at time 1 is uniform. In every
row the letters below the two tops are independent fair coins, independent of the tops. So the law at time 1 on
level 1 is π if and only if the pair (top of σ, top of τ) has the law: − with probability β for σ, fair for τ,
independent. Let a_(XY) be the total probability of the rows with (σ, τ) of type (X, Y). Then

    P(−, −) = a_FF/4 + a_(F,−F)/2,   P(−, +) = a_FF/4 + a_(F,+F)/2,
    a_(F,−F) = (1 − β)z²m₀/p + βz,   a_(F,+F) = (1 − β)z,   a_FF = (1 − β)z²/p + βz/m₀,
    a_(+F,−F) = p·a_(F,−F),   a_(+F,+F) = p·a_(F,+F).

So π is invariant if and only if

    (A)  (1 − β)z²m₀/p + βz = (1 − β)z            (then P(−,−) = P(−,+) and P(+,−) = P(+,+)),
    (B)  a_FF/4 + (1 − β)z/2 = β/2                (then P(−,±) = β/2, and P(+,±) = (1 − β)/2 because the total is 1).

(A): divide by z and use 1 − β = (1 + m₀)/(1 + 2m₀), 1 − 2β = 1/(1 + 2m₀): it reads z·m₀(1 + m₀)/p = 1, which is the
value of z.
(B): multiply by 4(1 + 2m₀)/z: (1 + m₀)z/p + 1 + 2(1 + m₀) = 2m₀/z. By (A) the first term is 1/m₀, and
2m₀/z = 4m₀³(1 + m₀) = 4m₀³ + 2m₀ + 1. So (B) reads 1/m₀ + 2 = 4m₀³, that is 4m₀⁴ = 2m₀ + 1. ∎

So the reduction of 625 to 32 cylinders is complete (the argument of its Step 1 is right), and in fact four
numbers per level suffice. `exact_check.py` verifies the invariance in Q(m₀) on the cylinders of depth 1, 2, 3, 4
(8, 32, 128, 512 cylinders), verifies (A) and (B), and shows that the check has teeth: with the active top fair, or
− with probability 1/3, the law is not invariant.

**4.2 The two martingales by hand.** c := e₁/(2(1 + m₀)) = 1/(1 + 2m₀), e₁ = 2(1 + m₀)/(1 + 2m₀), e₁m₀/(1 + m₀) = 2m₀/(1 + 2m₀).
ΔL¹ = ΔN_(σ+) + ΔN_(τ+) + e₁·ΔN_(τ−) − c·ΔΔ, ΔL² = ΔN_(σ−) − ΔN_(τ−) + ΔΔ/(1 + m₀), where ΔΔ = +1 for a go from level 1,
−1 for a go from level 0, 0 for a stay. From the table of 625 section 0: given a stay the active stack loses its
top + with probability 1/(1 + p) and, independently, the passive stack receives − with probability m₀/(1 + m₀); given a
go the passive stack receives + with probability p/(1 + p) and, if the active top is −, the active stack
independently loses it with probability 1/(1 + m₀). P(stay | top +) = 1/(1 + m₀).

| state | E[ΔL¹] | E[ΔL²] |
|---|---|---|
| level 1, top + | (1/(1+m₀))·[−2m₀/(1+2m₀) + 2m₀/(1+2m₀)] + (m₀/(1+m₀))·[1/(1+2m₀) − c] = 0 | (1/(1+m₀))·[−m₀/(1+m₀)] + (m₀/(1+m₀))·[1/(1+m₀)] = 0 |
| level 1, top − | 1/(1+2m₀) − c = 0 | −1/(1+m₀) + 1/(1+m₀) = 0 |
| level 0, top + | (1/(1+m₀))·[−2m₀/(1+2m₀)] + (m₀/(1+m₀))·[1/(1+2m₀) + c] = 0 | (1/(1+m₀))·[m₀/(1+m₀)] − (m₀/(1+m₀))·[1/(1+m₀)] = 0 |
| level 0, top − | −e₁/(1+m₀) + 1/(1+2m₀) + c = 0 | 1/(1+m₀) − 1/(1+m₀) = 0 |

The increments are bounded by B = 3. (`exact_check.py` verifies the four cases in Q(m₀).) The conceptual proof in
625 (derivative of the identity "the probabilities of the member u sum to 1" at the critical member, where the
derivative of log z vanishes) is also correct: ∂ log z/∂u₂ = 0 because p_− = q_−, and ∂ log z/∂u₁ = 0 is the equation
2p·m₀ = 1.

**4.3 Lemma 4 of 625.** Correct as written. For Theorem 2 only the following is used from Proposition 2: L¹ and L²
are martingales with increments bounded by B (items (b) and (c) are not needed: Doob's inequality with
E[(Lⁱ_n)²] ≤ nB² gives P(max_(k≤n) |Lⁱ_k| ≥ K′√n) ≤ B²/K′²). The other inputs are: Theorem 1 (the chain indexed by moves
is stationary under π, so at each fixed k the letters below the top are fair coins); Hoeffding's inequality with a
union bound over k ≤ n, j ≤ J_n and the two stacks, which gives P(G^c) ≤ 8(n + 1)·J_n·exp(−n^(0.8)/(2J_n)) → 0 for fixed K;
Lemma 3 (cancellation of letters pushed and popped), which is immediate; and (3.1), which I recomputed:
L¹ = x/2 + (1 + e₁)y/2 + D^σ/2 + (1 − e₁)D^τ/2 − cΔ, L² = x/2 − y/2 − D^σ/2 + D^τ/2 + Δ/(1 + m₀), det J = −(2 + e₁)/4 = −0.838.
Order of the constants: choose K′ with 2B²/K′² ≤ 1/4, put K = 2‖J⁻¹‖K′, then n₀ so large that P(G^c) ≤ 1/4 and
‖J⁻¹‖(√2·K′√n + C·n^(0.4) + C) < K√n for n ≥ n₀. Then P(max_(k≤n) max(|x_k|, |y_k|) ≤ K√n) ≥ 1/2 for n ≥ n₀.
Check (data, `mc_check.py`, 3,000 stationary runs of 2,000 moves): the frequency of the window event is 0.565,
0.916, 0.988, 1.000 for K = 1, 1.5, 2, 3; the level, the active top, the passive top and the second active letter
at time 2,000 have frequencies 0.505, 0.330 (β = 0.324), 0.506, 0.486 with standard error 0.009; Lemma 3 holds on
every run.

**4.4 Theorem 2 of 625, bookkeeping.**
- *Step (i), explicit.* From t every state (level, σ, τ) with bottom letter + in τ is reached in at most
  6(|σ| + |τ|) + 1 points, and t is reached from it in at most 2(|σ| + |τ|) + 2 points, by the moves named in 625;
  so C₁ = 6 works. I checked each named move against the table: the stay "σ receives −" at level 0 needs the top of
  τ to be +, which holds because τ = (+) while σ is built; the go "nothing changes" exists at every state; the stay
  "σ keeps, τ receives −" needs the top of σ to be +, which the procedure arranges; in the removal, a pop of the
  bottom letter of τ is never asked for. `macro_check.py` runs both procedures on all states with words of length
  at most 4 and on 300 random states with words up to length 40: every move is in the table.
  Every move of the true chain has probability at least δ = min(z², z²m₀, z²/p, z²m₀/p, z, zp, z/m₀, zp/m₀) > 0.
- *Step (ii).* Whether F occurs is decided by the level, the uniform variables and the top L letters of each stack:
  as long as the heights are at least −(L − 2), the active top is one of the top L − 1 initial letters or a pushed
  letter, and F fails at the first time a height leaves [−(L − 2), L − 2]. On F at least two letters of each window
  are never popped, so in the true chain started at s* the active stack always has a removable top (σ not empty
  at level 1; τ of length at least 2 at level 0), where by 3.2 the true kernel is the bulk kernel. So the true
  chain from s* performs the first T moves inside the window with probability at least 1/2.
- *Step (iii).* Let k* ≤ T be the first move index with T_(k*) ≥ T; T_(k*) ∈ {T, T + 1}. The set of paths counted is: the
  fixed path t → s*; a path γ of the true chain from s*, stopped at k*, inside the window; the fixed removal path
  from the end of γ to t; and pairs of gos "nothing changes" at t up to n points. Two different γ are not prefixes
  of each other (k* is a stopping time), so different γ give different paths from t to t with n points, and the
  probability of the union is the sum. The number of points left for the last part is even because every path
  t → t has an even number of points (stays have two points and the number of gos is even), and it lies between 0
  and 8C₁(K√n + 4) + 2 by the choice of T. So P_c(t → t in n points) ≥ δ^(C₁(2L+2)) · (1/2) · δ^(C₁(4L+4)) · z^(8C₁(K√n+4)+2)
  ≥ exp(−C√n) for all large even n. For the finitely many remaining m, M_m ≥ 1 (the path of gos "nothing changes"
  alone), so the constant C can be enlarged to cover every m ≥ 1.

## 5. Standard results used (none from the project)

The Jordan curve theorem for finitely many disjoint simple closed curves on the sphere (regions form a tree);
Doob's L² maximal inequality for martingales; Hoeffding's inequality for sums of independent fair signs.
No computer calculation is needed for the proof: every identity used is one of (E1) to (E4), (A), (B) or the table
of 4.2, all checked by hand above from 4m₀⁴ = 2m₀ + 1.
