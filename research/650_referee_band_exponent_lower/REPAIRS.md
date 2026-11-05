# 650 REPAIRS: arguments that are missing or thin in `643_band12_letter_balance/PROOF.md`, written out (2026-10-10)

Referee notes. Notation as in 643 `PROOF.md` section 0 and 632 `PROOF.md` section 0: the bulk chain with kernel P
on Ω = {0,1} × {+,−}^N × {+,−}^N (level, word σ, word τ, top letter first), its invariant product law π (625
Theorem 1), moves (index k) and points (T_k = number of points of the first k moves, a stay has 2 points, a go has
1), heights z_k = (x_k, y_k) (pushes minus pops), β = m₀/(1+2m₀), ω(+) = 1 − β, ω(−) = β, cyl(t) = {level 0, top of σ
is −, top of τ is +}, R(2n) = P_π(E_(2n) | cyl(t)). "Active" stack: σ at level 1, τ at level 0. Status of every item:
**proved** (complete argument in this file) unless marked otherwise. Scripts named here are in this folder.

Summary. I found no error that changes a statement. Sections 2 and 3 of 643 (the reversed machine, its
martingales, the pairing) are correct and I confirmed every exact claim with my own code. Theorem 4.1 is
correct; four places need an argument that 643 does not give (sections 1, 2, 4, 5 below), and one constant
changes (a factor 2 in the Claim of Step 3). No exponent changes.

## 1. The time reversal must be done on (state, move) cylinders

**What 643 does.** It defines P̂(y, x) = π(x)P(x, y)/π(y) as if Ω were countable, and proves Lemma 3.1 by
π(x_0)·ΠP(x_(i−1), x_i) = π(x_m)·ΠP̂(x_i, x_(i−1)). On Ω single states have π-measure 0, and a pair of states does
not determine the move (the word +++… is unchanged by a pop, for example), so both lines need a definition.

**Repair.** A *window* is w = (level; a_1…a_i; b_1…b_j) with i, j ≥ 1; [w] is the cylinder of states with that level
whose σ starts with a_1…a_i and whose τ starts with b_1…b_j. π[w] = (1/2)·ω(active top)·2^(−(i+j−1)).
A *forward move* is one of the ten entries of the move table (six at active top +, four at active top −); it
acts on windows (it reads the active top; it removes that letter or not; it puts a named letter on the passive
stack or not), and it maps [w] onto [mw] bijectively when w is long enough for the pop. Its probability
p(m | w) depends on (level, active top) of w. The *reversed move* m̂ undoes m: it removes the letter that m
pushed (so it needs that letter on top) and puts back the letter that m popped. The letter put back is
determined by the type of m: a stay pops +, a go pops −. So m̂ is a map on windows as well, [m̂(mw)] = [w].

*One-step identity.* For every window w of depth ≥ 2 in both stacks and every forward move m allowed at w:

    π[w]·p(m | w) = π[mw]·p̂(m̂ | mw),                                                   (1.1)

where p̂(m̂ | ·) is the table of Lemma R0 and depends on (level, top of σ, top of τ) of mw only. *Proof.* This is the
computation in the proof of Lemma R0 (ratio of two products of fair letters and one active top); I redid the
eight cases by hand and by machine (`ref_kernel.py` (1), exact in Q(m₀), all windows of depth 4: the ratio
π[w]p(m | w)/π[mw] is the same for all w with the same (level, tops) of mw and the same m̂, equals the table, and
the eight rows of the table sum to 1). ∎

*The reversed machine.* Define P̂ as the kernel of the stack machine that, at a state y, performs the reversed
move m̂ with probability p̂(m̂ | level, tops of y). By the row sums it is a probability kernel at every state of Ω
(not only almost every state). P̂_s is its law from s, P̂_π = ∫π(ds)P̂_s.

*Path identity.* For a window w_0 and forward moves m_1, …, m_k with w_i = m_i w_(i−1) (w_0 long enough that all
windows have depth ≥ 2), iterate (1.1):

    P_π(X_0 ∈ [w_0], the first k moves are m_1, …, m_k) = π[w_0]·Π p(m_i | w_(i−1)) = π[w_k]·Π p̂(m̂_i | w_i)
        = P̂_π(X̂_0 ∈ [w_k], the first k reversed moves are m̂_k, …, m̂_1).                                       (1.2)

On the left event X_k ∈ [w_k], and on the right event X̂_k ∈ [w_0]. Events of the form (1.2) generate the σ-field
of (X_0, m_1, …, m_k) on the left and of (X̂_0, m̂_k, …, m̂_1) on the right, and the map
(x_0, m_1, …, m_k) ↦ (x_k, m̂_k, …, m̂_1) is a bijection between the two path spaces. Hence:

**Lemma 3.1 (corrected statement).** For every k and every measurable set Γ of sequences (x_0, m_1, …, m_k)
(start state and k forward moves): P_π((X_0, m_1, …, m_k) ∈ Γ) = P̂_π((X̂_0 and the first k reversed moves) ∈ Γ̂),
where Γ̂ is the image of Γ under the bijection above. In particular (k = 0, 1) π is invariant for P̂ and P̂ is a version of the time
reversal of P with respect to π. (`ref_kernel.py` (1b) checks πP̂ = π directly on the 128 cylinders of depth 3.)

Heights, points, levels, "the floor is not popped", "the state after j moves lies in a set B" are functions of
(x_0, m_1, …, m_k), so every event used in 643 sections 3 and 4 is covered. The heights of the reversed path
relative to its start are the forward heights relative to the forward end.

**Which stack the reversed chain pops, and "fresh letters".** At level 1 (σ active): a reversed stay puts + on σ or
leaves σ alone (then it needs top of σ = +) and removes a − from τ or leaves τ alone; a reversed go (to level 0)
removes a + from σ or leaves it, and puts − on τ or leaves it. Level 0 is the mirror image. So the reversed
chain pops from both stacks at each level, a popped letter is + exactly when a reversed go pops it, and the law
of a reversed move depends on both tops. (`ref_kernel.py` asserts these facts for all 8 local states.) The
statement "a letter of the start state revealed for the first time is a fair coin" (625 Corollary 1) holds for
P̂_π with the same proof: the event "the first k reversed moves are a given sequence" is determined by the level,
the driving uniform variables and the start letters at depth ≤ (number popped so far) + 1 in each stack, and the
deeper letters are independent of these under π. It is not used anywhere in sections 3, 4 of 643.

**The local state is not a Markov chain, and this does not matter.** After a reversed pop the new top of that
stack is a hidden letter, so (level, top of σ, top of τ) is not Markov. Lemma R1 and Lemma R2 only need: the
conditional law, given the past, of the next move and of the next pair (level, active top) is a function of the
current triple (level, top σ, top τ). This holds: after a reversed stay the active top is + (kept, and then
required to be +, or pushed); after a reversed go the new active stack is the old passive one, which a reversed
go never pops, so its top is the old passive top or the pushed −. (`ref_kernel.py`, function `nxt`, asserts that
no hidden letter is needed.)

## 2. Lemma R3 written out: what 632 and 639 use, checked for P̂ item by item

643 proves Lemma R3 by "the same arguments". Here is the list of what the proofs of 632 Lemmas 1, 2 and of 639
sections 1, 3, 5.1, 5.2 use, and the check for the reversed chain. The reversed classes Â(r; ε, η, κ) and the
event Ĝ(r, ε, κ) are defined with P̂ and reversed move indices.

| Property used by 632 / 639 | Where used | For P̂ |
|---|---|---|
| (P1) π is invariant; the chain indexed by moves is stationary under π | Lemma 1 | section 1 above; `ref_kernel.py` (1b) |
| (P2) under π each word is one top letter and independent fair coins | Lemma 1 (Hoeffding) | same π |
| (P3) k ≤ T_k ≤ 2k | Lemma 1′ of 639; number of moves of a stage | reversed stay 2 points, reversed go 1 point |
| (P4) two martingales with increments ≤ B from every start state | Lemma 2 (b), (c) | Lemma R1; `ref_kernel.py` (2): drift 0 in all 8 local states, exactly; increments ≤ 1.36 |
| (P5) L = J·z + (1/2)δ^σD^σ + (1/2)δ^τD^τ + bounded, det J ≠ 0 | Lemma 2 (a) | (F3̂): a linear identity in the four net counts (checked, `ref_kernel.py`), det Ĵ = (2+e₁)/4 |
| (P6) D^σ_k = S(σ at time k; b_k) − S(σ at time 0; a_k) (625 Lemma 3) | Lemma 2 (a) | a statement about any process that pushes and pops letters on two stacks; holds verbatim |
| (P7) the predictable covariation is H·T + c·N + bounded + martingale with bounded increments | Lemma 2 (c) | Lemma R2; `ref_kernel.py` (3): the three functions lie in the span of the 7 drift vectors (rank 5), exactly |
| (P8) J⁻¹HJ⁻ᵀ = Σ | identification of the limit | Ĵ⁻¹HĴ⁻ᵀ = Σ exactly (`ref_kernel.py` (3)); also forced by reversal, see below |
| (P9) strong Markov property at the stopping times k(u) | stages | P̂ is a Markov kernel on Ω |
| (O1) Azuma, (O2) martingale functional central limit theorem | Lemma 2 (b), (c) | no property of the chain |

Not on the list, because 632 Lemmas 1, 2 and the one-stage bounds do not use them: the building and removal
procedures (F4), the identity (F1), the true chain. In 643 they enter for P̂ only through (3.1) (section 3 below).

**Lemma 1̂.** π(Â(r; ε, η, κ)^c) ≤ C_ε·(1 + κr²)·r·exp(−ε²r/4)·(1 + 1/η) for εr ≥ 4, with the constant of 632 Lemma 1.
*Proof.* π{a word is not ε-regular at scale r} ≤ C_ε·r·exp(−ε²r/4) by (P2), as in 632. By (P1),
P̂_π(Ĝ fails) ≤ 2(1 + κr²)·C_ε·r·exp(−ε²r/4), and π{s : P̂_s(Ĝ fails) > η} ≤ P̂_π(Ĝ fails)/η. ∎

**Lemma 1̂′.** P̂_π(F ∩ {X̂_(k̂(u)) ∉ A}) ≤ (u/2 + 1)·π(A^c) for every event F, every measurable A and integer u ≥ 0, k̂(u) the
first reversed move with T̂_k ≥ u. *Proof.* (P3) gives ⌈u/2⌉ ≤ k̂(u) ≤ u; union bound and (P1). ∎

**Lemma 2̂ (a).** Let s ∈ Â(r; ε, η, κ), ε ≤ ε₀. On Ĝ, for all k ≤ κr²: Ẑ*_k ≤ C·L̂*_k + C·εr + C and
abs(ẑ_k − Ĵ⁻¹L̂_k) ≤ C·ε·(r + L̂*_k) + C. *Proof.* By (P6) and regularity of the words at reversed times 0 and k,
abs(D̂^σ_k) ≤ ε·max(b̂_k, r) + ε·max(â_k, r) ≤ ε·(2r + 3Ẑ*_k), the same for τ. Insert in (P5) and solve, as in 632. ∎

**Lemma 2̂ (b).** P̂_s(max over k ≤ k̂(κr²) of abs(ẑ_k) ≥ λr) ≤ 4·exp(−cλ²/κ) + η for λ ≥ 1, r ≥ r₀. *Proof.* (a), (P3), (P4), (O1). ∎

**Lemma 2̂ (c).** For every κ, ε′ > 0 there are ε, η₀, r₀ such that for r ≥ r₀ and s ∈ Â(r; ε, η₀, κ) the law under P̂_s of
u ↦ r⁻¹·ẑ_(k̂(r²u)), 0 ≤ u ≤ κ, is within Prokhorov distance ε′ (uniform metric) of the law of W, the Brownian motion
with covariance Σ started at 0.
*Proof.* Suppose not: there are r_n → ∞, ε_n, η_n → 0, s_n ∈ Â(r_n; ε_n, η_n, κ) at distance > ε′. For integers v put
Ĝ_v := F̂_(k̂(v)). As in 639 section 3: k̂(v+1) − k̂(v) ∈ {0, 1}, and it is 1 exactly when T̂_(k̂(v)) = v, an event of Ĝ_v. So
v ↦ L̂_(k̂(v)) is a martingale for (Ĝ_v) by (P4) (optional stopping at bounded stopping times), with increments
≤ 3, and its predictable covariation is Â_(k̂(v)), the process of Lemma R2. By Lemma R2,
Â^(ij)_k − H_ij·T̂_k = c_ij·N̂(k) + (g′_ij(X̂_k) − g′_ij(X̂_0)) + R̂^(ij)_k. On Ĝ, abs(N̂(k)) ≤ abs(ẑ_k) + abs(D̂^σ_k) + abs(D̂^τ_k) ≤ C·(r_n + Ẑ*);
Ẑ*/r_n is tight by (b); max abs(R̂)/r_n is tight by Doob (R̂ has bounded increments and there are at most κr_n² moves);
T̂_(k̂(v)) ∈ {v, v+1}; P̂(Ĝ fails) ≤ η_n → 0. So r_n⁻²·Â_(k̂(r_n²u)) → H·u in probability for each u, the jumps of the
martingale are ≤ 3/r_n and those of the covariation ≤ 9/r_n². By (O2) in the form checked in 639 section 3,
r_n⁻¹·L̂_(k̂(r_n²·)) converges in law to a Brownian motion with covariance H, and by (a)
sup_u abs(r_n⁻¹ẑ_(k̂(r_n²u)) − Ĵ⁻¹·r_n⁻¹L̂_(k̂(r_n²u))) → 0 in probability. The limit of the heights is Ĵ⁻¹·(Brownian motion with
covariance H), which has covariance Ĵ⁻¹HĴ⁻ᵀ = Σ by (P8). Contradiction. ∎

Remarks. (i) (P8) is also forced without computation: under P_π and P̂_π the height paths are time reversals of
each other, so they have the same covariance at every time; the exact computation confirms it, and it is
needed as stated because Lemma 2̂ (c) is about fixed start states, where H enters through Lemma R2.
(ii) Ĵ is J with its columns exchanged; since Σ₁₁ = Σ₂₂, (P8) for Ĵ is the same identity as for J.
(iii) `ref_mc.py` (data) runs P̂ from one fixed non-random start state and from π; the covariance per point of the
heights agrees with Σ in both cases (section 7).

## 3. The pairing (Lemma 3.1, (3.1), Proposition 3.2, Corollary 3.3): checked, with the points of the brief

**(a) True chain against bulk chain; how t is embedded.** R(2n) = P_π(E_(2n) | cyl(t)) is 625 (2.1), refereed in 630
(3.3 and 4): on cyl(t) ∩ {both heights ≥ 0} the bulk chain is the true chain started at t; the top − of σ stands for
"σ empty", the top + of τ for the bottom letter of τ, and a path dies exactly when one of the two is popped. The
kernel of the true chain is already the transformed one (weights z^(points)·V(s′)/V(s)); a path t → t with 2n
points has weight z^(2n), and no further eigenfunction factor appears. My program reproduces the integers
M_m = Λ^(m−1)·R(2m − 2) = 1, 2, 7, 31, 161, 932, 5841, 38891 exactly in Q(m₀) (`ref_pairing.py` (a)), so the table,
the boundary rule and the identity are the right ones.

**(b) Point-time against move index.** The split is at a move boundary: the paths with T_k = n for some k. The
paths in which a stay covers the points n, n + 1 are dropped (a lower bound). Lemma 3.1 is applied to each k
separately, and the events {T̂_k = N} are disjoint for different k, so (3.1) is an equality. Stationarity is used
in the move index only: (1.2) is an identity between a forward path of k moves and a reversed path of k moves.

**(c) The identity (3.1) and Proposition 3.2.** Let (s, z) be a bulk state with a floor (letters −, + at depths
z_x + 1, z_y + 1) and C(s, z) its cylinder (level and the top z_x + 1, z_y + 1 letters). The left side of (3.1) is
π(cyl t)·G_N(s, z), G_N(s, z) the probability that the true chain started at t is, after exactly N points, at the
true state "top z letters of s above t". Under the bijection of section 1 the forward paths of this event
correspond to the reversed paths that start in C(s, z), never pop a floor letter and end, after exactly N
points, at level 0 with both heights −z. So (3.1) holds: π(cyl t)·G_N(s, z) = π(C(s, z))·ê_N(s, z). Markov property at k:

    R(2n)·π(cyl t) = Σ_(true states (s,z)) π(cyl t)·G_n(s, z)·e_n(s, z) + (straddling paths) ≥ Σ_((s,z)) π(C(s, z))·ê_n(s, z)·e_n(s, z).

For fixed z the cylinders C(s, z) partition {s : z is a floor of s}, and ê_n, e_n are constant on each, so the right
side is Σ_z ∫π(ds)·1[z is a floor of s]·ê_n(s, z)·e_n(s, z). Direction and measurability are in order.
*Exact test* (`ref_pairing.py`, Q(m₀), my forward table and the reversed kernel computed from π, N ≤ 7, up to 120
states): π(cyl t)·G_N(s) = π(C(s))·ê_N(s) for every reachable true state, and R(2N) = Σ_s G_N(s)e_N(s) + straddling
part, both as exact identities. The pairing sum is 50 to 86 percent of R(2N) for N ≤ 7.

**(d) Jensen.** For z ∈ box(√n) the floor letters are at depth ≥ 2, hence fair and independent of the level:
π{z is a floor of s, level ≡ n} = 1/8. So the part of the sum with z ∈ box(√n) equals (#box/8)·E_(U_n)[ê_n·e_n], and
E[exp(log ê_n + log e_n)] ≥ exp(E log ê_n + E log e_n). Both logarithms are bounded on the support of U_n for
n ≥ n_f (a constant): by (4.1), −log e_n and −log ê_n are at most C_f·√n + B·log n there. So the means are finite and the
inequality is not empty. The factor n in R ≥ c·n·exp(…) is #box(√n) = (3√n + O(1))².

## 4. The Claim of Step 3: reversed stopping times, and a factor 2

**The point.** κ_k is a stopping time of the forward chain; read backward, "κ_k = m" is not a stopping event, and
one reversed trajectory can be the image of two forward paths of different lengths.

**Corrected Claim.** For k_θ < k ≤ J and every measurable set B of bulk states,

    Σ_(z∈box(√n)) P_π(F_1(z) ∩ … ∩ F_k(z) ∩ {X_(κ_k) ∈ B}) ≤ 2·#box(r_k)·[ π(B)·q_S^(k_θ) + k_θ·(n+1)·max_(j≤k_θ) π(Â^g_j c) ].

*Proof.* Let D := n − ⌊r_k²⌋. {κ_k = m} = {T_(m−1) < D ≤ T_m}. Under the bijection of section 1 (forward path of m moves
from X_0 to X_m ↔ reversed path of m moves from X̂_0 = X_m), T_m = T̂_m and T_(m−1) = T̂_m − (points of the first
reversed move). So the image of {κ_k = m} is the event Q_m := {T̂_m = D} ∪ {T̂_m = D + 1 and the first reversed move
is a stay}; I say that m *qualifies*. For a floor z put z′ := z + (heights at κ_k); on F_k(z), z′ ∈ box(r_k), and z ↦ z′
is injective. The image of F_1(z) ∩ … ∩ F_k(z) ∩ {κ_k = m} ∩ {X_(κ_k) ∈ B} is contained in

    {X̂_0 ∈ B} ∩ Q_m ∩ Γ̂_m(z′),   Γ̂_m(z′) := {reversed heights ≥ −z′ up to move m; for each j < k the floor
                                    z′ + (reversed heights) lies in box(r_j) at the reversed image m_j of κ_j},

where m_j is the reversed move index with T̂_(m_j) = T_(κ_k) − T_(κ_j) ∈ [D_j − 1, D_j + 1], D_j := ⌊r_j²⌋ − ⌊r_k²⌋ (both
overshoots are 0 or 1). By Lemma 3.1 (section 1), summing over m and z,

    left side ≤ Σ_(z′∈box(r_k)) Ê_π[ 1(X̂_0 ∈ B)·#{m : Q_m ∩ Γ̂_m(z′) holds} ].

At most two m qualify (one with T̂_m = D, one with T̂_m = D + 1). `ref_stop.py` confirms the identity behind this
line exactly for D = 5, 7 and shows that two do occur: the sum over m is 1.17 times the probability of the union
in those examples. Now define the stopping times of the reversed chain ς_j := first reversed move with at
least D_j − 1 points, 0 ≤ j ≤ k_θ. Then T̂_(ς_j) ∈ {D_j − 1, D_j}, so ς_j ≤ m_j ≤ ς_j + 2 and, on Γ̂_m(z′), the floor at ς_j lies in
box⁺(r_j) = [r_j − 2, 4r_j + 2]². Also ς_0 ≤ m for a qualifying m (T̂_m ≥ D = D_0), and between ς_j and ς_(j−1) there are
at least D_(j−1) − 1 − D_j ≥ (1/δ − 1)·r_j² − 3 points (r_(j−1)² = r_j²/δ for j ≤ k_θ). Hence, for every z′ and m,

    {X̂_0 ∈ B} ∩ Q_m ∩ Γ̂_m(z′) ⊂ E(z′) := {X̂_0 ∈ B} ∩ ⋂_(j=1..k_θ) {floor at ς_j in box⁺(r_j); the floor z′ is not popped between ς_j and ς_(j−1)},

and E(z′) does not depend on m. So the left side is at most 2·Σ_(z′) P̂_π(E(z′)). On the event that X̂_(ς_j) ∈ Â^g_j for all
1 ≤ j ≤ k_θ, the strong Markov property of P̂ at ς_1, then ς_2, …, ς_(k_θ) (the indicator of {X̂_0 ∈ B} is measurable at
time 0) and (St3) for P̂ at each of them give P̂_π(E(z′) ∩ all regular) ≤ π(B)·q_S^(k_θ). The complement has probability
at most Σ_j (D_j/2 + 1)·π(Â^g_j c) ≤ k_θ·(n + 1)·max_j π(Â^g_j c) by Lemma 1̂′ (ς_j = k̂(D_j − 1)). ∎

The factor 2 multiplies E_(U_n)P̃(ϑ = k) in Step 4 (i) and disappears in C_δ. Nothing else changes.

## 5. Theorem 4.1: the remaining small points

**5.1 (4.1), ranges.** After stage k the number of points to go is N_k = ⌊r_k²⌋ − e, e ∈ {0, 1} the overshoot, so
r_k² − 2 < N_k ≤ r_k²; (4.1) is stated for r² − 1 ≤ N ≤ r² + 1. Its proof gives the same bound for r² − 2 ≤ N ≤ r² + 1
(only N ≥ L, N ≥ L′ and the parity are used: L ≤ 16r + 12, L′ ≤ 48r + 31, so r ≥ 50 suffices). Parity: a point changes
the parity of the level exactly when it belongs to a go, so level(X_(κ_k)) ≡ level(s) + (n − N_k) ≡ N_k mod 2 when
level(s) ≡ n; the floor of the start is a floor of X_(κ_k) on F_1 ∩ … ∩ F_k (floor letters are not popped). So (4.1)
applies at every stage end. For ê the bound has the extra term log(1/π(cyl t)) (a constant; enlarge C_f).

**5.2 The chain P̃.** Its state is (bulk state, floor, overshoot e ∈ {0, 1}); the k-th step from (x, z, e) is the law
under P_x of (X_κ, z + heights at κ, overshoot) given F_k, κ the first move with at least ⌊r_(k−1)²⌋ − ⌊r_k²⌋ − e points.
It is only needed at states with x ∈ A_(k−1) and z ∈ box(r_(k−1)), where P_x(F_k) ≥ q_k > 0 by (St1), (St2) (these allow
abs(D − (1 − δ)r²) ≤ 3). The bound P̃_((s,z))(ϑ = k) ≤ (q_1⋯q_k)^(−1)·P_s(F_1 ∩ … ∩ F_k ∩ {X_(κ_k) ∉ A_k}) is correct: on
{ϑ ≥ j} the j-th conditioning divides by P_(X_(κ_(j−1)))(F_j) ≥ q_j.

**5.3 Where the logarithm meets a probability.** Only in Step 1: −log E_s[e_(N_k) | F_k] ≤ E_s[−log e_(N_k) | F_k] (Jensen for
the conditional law, legitimate because P_s(F_k) > 0 and e_(N_k) > 0 on F_k by (4.1)). After that every step is
linear in probabilities: V_0 ≤ Σ log(1/q_k) + β(r_J) + Σ_k β(r_k)·P̃(ϑ = k) is a pointwise inequality in (s, z), and
E_(U_n) is applied to both sides. No expectation of a logarithm is exchanged with a conditional probability.

**5.4 No circularity.** The only lower bound for return probabilities that enters is Theorem Lo of 632
(R(2j) ≥ j^(−B), all j ≥ 1; refereed in 639), through the fallback (4.1). The bound on the conditioned chain uses
upper bounds only ((St3), Lemma 1, Lemma 1′ and their reversed forms).

**5.5 Order and independence of the constants** (checked at the sources). Absolute: δ_* = 1/16; q_*, ε_*, η_*, r_* (from
Lemma 2 (c) and 2̂ (c) with κ = 1 and one fixed Brownian event; take the worse of the two chains); B (Theorem Lo);
C_f (from δ₀, B, π(cyl t)); c₄ ((BM4): a Brownian constant of the box [1, 4]², see section 6). Then: ε; θ and ε₁ with
(1 + 32C_f/ε_*²)·2θ·log(1/q_*)/log 16 + 16C_f·ε₁/ε_*² ≤ ε/4; C₁ = C₁(ε₁) (639 Proposition B1 with μ′ = μ − ε₁/2; it does
not depend on δ: q_S ≤ C·((6 + δ)²·δ/(1 − 2δ))^(μ′) + δ^(1+μ)); δ with (1 + 32C_f/ε_*²)·log(2C₁/c₄)/log(1/δ) ≤ ε/4; then
(ε_r, η₀, r₀) and (ε_g, η_g, r_g) for both chains; then n₀. With these the displayed condition of 643 holds.
ϰ ≥ ε₁/2 > 0, so r^♯ → ∞ and all thresholds r_*, r₀, r_g, C_f are passed for n ≥ n₀.

**5.6 The exponent.** Σ_k log(1/q_k) ≤ k_θ·((1+μ)·log(1/δ) + log(2/c₄)) + (J − k_θ)·log(1/q_*), k_θ ≤ log n/log(1/δ),
J − k_θ ≤ (2θ·log n + log(1/δ))/log 16: this is (1 + μ + o(1))·log n per half. With #box(√n) ≥ n the two halves give
R(2n) ≥ c·n·n^(−2(1+μ+ε)) = c·n^(−α−2ε), α = 1 + 2μ. The count is right: each half pays μ for survival in the wedge and
1 for hitting the apex (area δ per stage), and the factor n is the number of floors.

## 6. (BM4) of 638 (the one unrefereed outside input): checked

Statement used in (St1): for a = 1, b = 4 there is c₄ > 0 such that for all δ ∈ (0, 1/8), ζ ∈ [1, 4]², θ′ ≤ √δ/8 and all
u ∈ [3/4, 1]: P_ζ(W ∈ Q + θ′(1,1) on [0, u], W_u ∈ [√δ + θ′, 4√δ − θ′]²) ≥ c₄·δ^(1+μ). The proof in 638 section 2 is correct:
after the translation by θ′(1,1) the start lies in [1/2, 4]², the target contains √δ·D₁ with D₁ = [5/4, 3]²; by
symmetry of the killed kernel p^Q (in the coordinates in which W is standard) and Chapman–Kolmogorov,
p^Q_u(ζ′, y) ≥ inf_(w∈D′) p^Q_(u/2)(w, ζ′)·P_y(τ_Q > u/2, W_(u/2) ∈ D′); the infimum is over a compact set of (w, ζ′, u) in
the open wedge, positive and free of δ; the last factor is ≥ c·abs(y)^(2μ) ≥ c·δ^μ by 639 Proposition B2 (Z = D₁,
λ = √δ, time u/2 ∈ [3/8, 1/2], which is covered after scaling time by 2); the area of √δ·D₁ is 3.0625·δ. So c₄ does not
depend on δ or u. (BM4) inherits the status of B2: proved from the classical heat-kernel series of the wedge
(639 (2.1), status literature, tested numerically there). With Lemma 2 (c), accuracy min(θ′, c₄δ^(1+μ)/2), this
gives (St1) for P, and with Lemma 2̂ (c) for P̂.

## 7. Data (not part of any proof)

`ref_mc.py`: the reversed machine (kernel computed from π and the forward table) run for U points, from π and
from the fixed state (level 1, σ = +−+−…, τ = −+−+…). Var x/U, Var y/U, Cov/U: 0.347, 0.337, −0.146 (π) and 0.344,
0.345, −0.136 (fixed) at U = 3000, 4000 runs; 0.343, 0.364, −0.144 (π) and 0.364, 0.340, −0.147 (fixed) at
U = 30000, 3000 runs; Σ = (0.3502, 0.3502, −0.1489), standard error about 0.009. Consistent with Lemma 2̂ (c).
