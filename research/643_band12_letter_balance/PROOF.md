# 643: band [−1, 2]. The lower bound M_m ≥ Λ^m·m^(−α−o(1)) without (DR*): time reversal, a pairing at half time, and Jensen (2026-10-10)

Status words as in `WAVE_620_COMMON.md`. **(DR*) of 638 is NOT proved here** (section 1 says what was found and why
the plan "balanced revealed letters give balanced pushed letters, by Azuma" does not work as posed). **What is proved
here (complete arguments, not reviewed) is the theorem that (DR*) was meant to give, by a route that does not use
(DR*), Proposition G or Lemma 2′ of 638:**

    M_m = Λ^m·m^(−α+o(1)),   α = 1 + π/arccos((1+2m₀)/(3+4m₀)) = 3.7757713065…      (Theorem 5.1)

The new inputs are: (i) the time reversal P̂ of the bulk chain with respect to π is an explicit stack machine with
the same two martingales, the same quadratic-variation identity and the same covariance (section 2; exact
arithmetic in Q(m₀), `rev_exact.py`); (ii) an exact pairing at half time, R(2n)·π(cyl t) ≥ Σ π·ê_n·e_n, in which
both factors are descent probabilities from the same state, one for P and one for P̂ (section 3); (iii) Jensen's
inequality, which reduces the pairing to the *typical* descent probability from a stationary state, for each of the
two chains (section 3); (iv) for a typical start the conditioned chain stays typical with respect to π up to a
factor n^(o(1)), because the density of its law is a growth probability of the *other* chain (section 4).

## 0. Notation and inputs

All notation of 632 `PROOF.md` section 0 and 638 `PROOF.md` section 0: m₀, Λ, the bulk chain with kernel P on
Ω = {0,1} × {+,−}^N × {+,−}^N, π (625 Theorem 1), P_s, P_π, moves and points, T_k, k(u), heights z_k = (x_k, y_k),
cyl(t) = {level 0, top of σ is −, top of τ is +}, Q the quarter plane, W the Brownian motion with covariance Σ,
μ = π/(2θ) = 1.3878856532, α = 1 + 2μ, δ₀ the smallest move probability, the classes A(r; ε, η, κ), (F1) to (F4),
R(2n) = P_π(E_(2n) | cyl(t)), M_m = Λ^(m−1)·R(2m − 2).
Constants of the move table: a := 1/(1+m₀) (probability of a symbol of class A), w := 2m₀/(1+2m₀) (coin "lose"),
λ₀ := m₀/(1+m₀) (coin "receive −"), p₊ := 1/(1+2m₀) (coin "receive +"), β := m₀/(1+2m₀), ω(+) := 1 − β, ω(−) := β
(the law of the active top under π).

*Floors.* For a bulk state s and z = (z_x, z_y) ∈ N², z is a *floor* of s if the letter of σ at depth z_x + 1 is −
and the letter of τ at depth z_y + 1 is +. For a floor z and N ≥ 0,

    e_N(s, z) := P_s(heights ≥ −z during N points; for some k: T_k = N, level 0 and heights = −z after k moves),

and ê_N(s, z) is the same probability for the reversed chain P̂ of section 2. Both depend on s only through its
level and its top z_x + 1 and z_y + 1 letters; C(s, z) denotes the cylinder of the states that agree with s there.
e_N(s, z) is the probability that the true chain started at the state "the top z letters of s above t" is at t
after exactly N points. Parity: e_N(s, z) > 0 requires N ≡ level(s) mod 2; the same for ê.

Inputs, with their status. From 632 as refereed in 639: Lemma 0, Lemma 1, Lemma 1′ of 639, Lemma 2 (a), (b), (c),
the stage bound "one stage, near start" of 632 section 2.1 with Proposition B1 of 639, Theorem Lo
(R(2n) ≥ n^(−B)), Theorem U, (F1) to (F4), (O1), (O2). From 638 (not refereed): only (BM4) (section 2 there; it rests
on Proposition B2 of 639 and on positivity and continuity of the killed heat kernel of a wedge). 625 Theorem 1.
Not used: (DR*), Proposition G, Proposition T, Proposition D, Lemma 2′, Theorem A, Theorem B of 638.

## 1. The pushed layer: exact description, the linear response, and why (DR*) was not closed

### 1.1 Rounds (the autonomous description asked for in task 1)

A *visit* is a maximal run of moves at one level; it consists of stays and ends with a go. During a visit at
level 1 the stack σ only loses letters and τ only receives. With the driving symbols of 638 section 5.1, a visit at
level 1 with σ = s₁s₂… (top first) is described as follows. While the top is +: with probability
g := (1 − a)/(1 − a + a·w) = (1+2m₀)/(3+2m₀) = 0.58636 the visit ends at this letter (a go at top +, the letter stays);
with probability 1 − g = 2/(3+2m₀) the letter is removed and the next letter becomes the top. If the top is −, the
visit ends there (a go at top −), and the letter is removed with probability a. So **a visit eats a geometric run
of leading + letters and stops at the first −.** The number ν of stays of the visit is the number of removed +
letters plus an independent geometric number of idle stays for each + letter that was on top. The visit writes on
τ the block (−)^n (+)^e, n binomial (ν, λ₀), e Bernoulli (p₊) independent of everything else. The new top of τ is
tested at once by the next visit (level 0), which is the erase phase of τ and the write phase of σ.

*Answer to "which finite information".* None is finite. Given the i.i.d. symbols, the pushed letters are a
function of the sequence of the letters that are on top at the symbols of class A, and these are pushed letters
or letters of the start state. The two stacks interact only through the numbers ν of stays of the visits.

### 1.2 Letter fates and the exact linear response

Let R_(σ±), W_(σ±) be the gross numbers of letters ± removed from and written on σ during the first k moves, the same
for τ, G^(1±), G^(0±) the go counts of 638 Lemma F, and St¹, St⁰ the numbers of stays at level 1, 0. From the move
table, each of the following is a martingale with bounded increments, from every start state:

    G^(1+) − ((1+2m₀)/2)·R_(σ+),   G^(1−) − (1+m₀)·R_(σ−),   St¹ − G^(1+)/m₀,   W_(τ−) − λ₀·St¹,   W_(τ+) − p₊·(G^(1+) + G^(1−)),

and the mirror images. (Per move at level 1 and top +: go with probability 1 − a, removal with probability a·w,
ratio (1+2m₀)/2; at top −: one go, removal with probability a = 1/(1+m₀); stays and gos at top + have ratio 1/m₀.)
Eliminating the go counts gives the **gross flow identity**: with κ := (1+m₀)/(1+2m₀) = 0.67636,

    (W_(τ+), W_(τ−))ᵀ = 𝐌·(R_(σ+), R_(σ−))ᵀ + martingale,    𝐌 = [[1/2, κ], [1/(2κ), 0]],    eigenvalues 1 and −1/2.

Left eigenvectors: (1, κ) for 1 and (1, −2κ) for −1/2. So Φ := #(+) + κ·#(−) is conserved by the transfer from one
stack to the other (this is the martingale L¹ + κL² of 625), and **Ψ := #(+) − 2κ·#(−), the deviation of a flow
from the equilibrium composition (2κ : 1), is multiplied by −1/2 at each transfer: the linear response has gain
1/2 and a sign change.** This is the exact form of the "bounded gain" asked for in task 1. Under π the gross flows have
57.5 percent of letters + (2κ/(1+2κ)), while the letters that survive are fair. These identities are equivalent to
Lemma F of 638 (same four relations); they say nothing about which letters survive.

### 1.3 Why task 2 does not work as posed (three exact reasons), and what was not proved

(a) *The response is not a function of window sums.* The law of the block written during a visit depends on the
run structure of the letters on top, not on their density: P(the visit ends at a +) is g/(1+g) = 0.3696 for a
visit that starts on fresh fair coins, and in the stationary regime of visits on the pattern (+−)^∞ with nothing
written on top it is π₊·g = 0.4250 (π₊ = 1/(1 + m₀(1−g)) = 0.7249). Both inputs are perfectly balanced in windows.
(b) *The martingale errors of Lemma F are of the order of the window.* The top b letters of σ at time k have sum
D^σ_k − D^σ_t, t the last time at which the height was x_k − b; typically k − t is of order b², so the martingale
terms of Lemma F are of order b, not εb. Balance of a window is finer than any identity between go counts.
(c) *The exact likelihood ratio.* A removed letter is + if it is removed by a stay and − if by a go, so a move
sequence determines the letters it removes. Hence P_s(γ) = P_π(γ)/π(C_γ(s)) for every move sequence γ compatible
with s, where C_γ(s) fixes the revealed letters: Proposition D of 638 is an identity, not only a bound, and
P_s = P_π(· | the letters of the start state). So (DR*) asks for a bound that is uniform over *every* regular start
word, for a conditional law of the stationary chain. I found no argument that uses the revealed letters through
their statistics only, and no counterexample. Heuristically (DR*) is true: of the letters of σ that survive, only
O(1) per dyadic scale are written at moves at which a letter of the start state of τ is on top.
**Status of (DR*): conjecture, unchanged.** The observation that led to the route below: a first irregular window
must contain the letter just pushed, and "for every regular start" can be replaced by "for a π-typical start in
logarithmic mean" as soon as the two halves of the excursion are paired by Jensen's inequality (section 3).

## 2. The time reversal of the bulk chain

P̂(y, x) := π(x)·P(x, y)/π(y) is the kernel of the stationary bulk chain read backward; π is invariant for P̂.
P̂_s, P̂_π are its laws. A *reversed move* has the points of the forward move it reverses (stay 2, go 1). For a state
y at level 1 let A be the top of σ and P the top of τ (at level 0 exchange the stacks; "active" and "passive" below
always mean the roles in y for the forward chain).

**Lemma R0 (the reversed machine; proved, exact).** The reversed moves from y and their probabilities are

    reversed stay (level kept, 2 points; the forward predecessor had active top +):
      (S1) nothing changes                                  needs A = +            a(1−w)(1−λ₀)
      (S2) the passive stack loses its top −                needs A = +, P = −     2a(1−w)λ₀
      (S3) the active stack receives +                                             (ω(+)/(2ω(A)))·a·w·(1−λ₀)
      (S4) active receives +, passive loses its top −       needs P = −            (ω(+)/ω(A))·a·w·λ₀
    reversed go (level changes, 1 point):
      (G1) nothing changes                                                         (ω(P)/ω(A))·(1−a)(1−p₊)
      (G2) the active stack loses its top +                 needs A = +            (2ω(P)/ω(+))·(1−a)·p₊
      (G3) the passive stack receives −                                            (β/(2ω(A)))·a·(1−p₊)
      (G4) active loses its top +, passive receives −       needs A = +            (β/ω(+))·a·p₊

The law of a reversed move depends on (level, A, P) only. After a reversed stay the active top is +; after a
reversed go the new active stack is the old passive one, and its top is P or the pushed −. So the law of the next
(level, active top) and of all increments used below is a function of the *local state* (level, s₁, t₁): 8 states.

*Proof.* π(x)/π(y) is computed from the product form of π: a fair letter gained by x counts 1/2, and the active
tops count ω. The eight cases are the forward moves "stay" (four) and "go" (four, the two tops of the predecessor
give the same probability (1−a) for a go without pop). `rev_exact.py` (1) computes π(x)P(x, y)/π(y) on windows of
depth 4 from the forward table in Q(m₀), finds that it depends on (level, s₁, t₁) only, that the eight rows sum to
1, and that it equals the table above. ∎

Remark. P̂ is not the forward machine after exchanging the stacks and flipping the letters (638 section 5.2): the
conditions on the tops hold after the move instead of before, and both tops enter the probabilities. It is a
machine of the same kind: a letter revealed by a pop for the first time is fair (625 Corollary 1 holds for P̂ with
the same proof), a popped letter is + if popped by a reversed go and − if by a reversed stay.

For the reversed chain let N̂_c(k), c ∈ {σ+, σ−, τ+, τ−}, be the net letter counts (pushed minus popped during the first
k reversed moves), x̂ = N̂_(σ+) + N̂_(σ−), ŷ = N̂_(τ+) + N̂_(τ−) the heights, T̂_k the points, e₁ = 2(1+m₀)/(1+2m₀).

**Lemma R1 (two martingales; proved, exact).** Let g₁ := −p₊, +p₊ at (level 1, active top +), (level 1, active top −)
and 0 at level 0; g₂ := a at (level 0, active top −) and at (level 1, active top +), 0 otherwise. Then

    L̂¹_k := e₁·N̂_(σ+) + N̂_(σ−) + N̂_(τ−) + g₁(X̂_k) − g₁(X̂_0),        L̂²_k := N̂_(τ+) − N̂_(σ+) + g₂(X̂_k) − g₂(X̂_0)

are martingales of the reversed chain from every start state, with increments bounded by 3.
(`rev_exact.py` (2): the drift is 0 in the 8 local states, exactly. The letter parts are the mirror images of
L¹, L² of 625 under exchange of the stacks and flip of the letters.)

**(F3̂)** With D̂^σ_k := N̂_(σ+) − N̂_(σ−) = S(σ at reversed time k; b̂_k) − S(σ at reversed time 0; â_k) (625 Lemma 3; it is
a statement about any process of two stacks) and the same for τ:

    L̂_k = Ĵ·ẑ_k + (1/2)·((e₁−1), −1)ᵀ·D̂^σ_k + (1/2)·(−1, 1)ᵀ·D̂^τ_k + (a term bounded by 2),   Ĵ = [[(1+e₁)/2, 1/2], [−1/2, 1/2]],   det Ĵ = (2+e₁)/4.

**Lemma R2 (Lemma 0 of 632 for the reversed chain; proved, exact).** With the same matrix H as in 632 Lemma 0 there are constants c_ij ∈ R⁴ and functions g′_ij of (level, active top) such that
Â^(ij)_k := Σ_(l<k) E[ΔL̂ⁱ_(l+1)·ΔL̂ʲ_(l+1) | F̂_l] = H_ij·T̂_k + c_ij·N̂(k) + g′_ij(X̂_k) − g′_ij(X̂_0) + R̂^(ij)_k, R̂ a martingale with
bounded increments. Moreover Ĵ⁻¹·H·Ĵ⁻ᵀ = Σ, the covariance of 625 Theorem 3.

*Proof.* `rev_exact.py` (3): the functions f = E[ΔL̂ⁱΔL̂ʲ | ·] − h·E[ΔT̂ | ·] of the local state, h = ω̂·f/ω̂·E[ΔT̂ | ·] (ω̂ the
π-weights of the 8 local states), are exact linear combinations of the four drift vectors E[ΔN̂_c | ·] and of the
three coboundaries P̂g − g, g an indicator of a pair (level, active top). (The seven vectors have rank 5; the three
right sides lie in their span, which the program checks by exact elimination.) Sum along the path. The values
h are H₁₁, H₁₂, H₂₂ of 625, and Ĵ⁻¹HĴ⁻ᵀ has the entries of Σ, exactly. ∎

**Reversed classes.** Â(r; ε, η, κ) := {s : both words of s are ε-regular at scale r and P̂_s(at some reversed move
k ≤ κr² a word is not ε-regular at scale r) ≤ η}.

**Lemma R3 (Lemmas 1, 1′, 2 of 632 and 639 for the reversed chain; proved by the same arguments).**
(1̂) π(Â(r; ε, η, κ)^c) ≤ C_ε·(1 + κr²)·r·exp(−ε²r/4)·(1 + 1/η) for εr ≥ 4; P̂_π(F ∩ {X̂_k ∉ Â}) ≤ π(Â^c) at every fixed
reversed move index k, and with the factor (u/2 + 1) at a point-time u.
(2̂) Lemma 2 (a), (b), (c) of 632 hold for P̂, the classes Â and the same limit W (Brownian motion with covariance Σ).

*Proof.* The proofs of 632 Lemmas 1, 2 and of 639 sections 1, 3, 5.2 use only: π is invariant and is a product of
a top letter and fair coins (Lemma 1); two martingales with bounded increments from every state (here Lemma
R1); the identity (F3) with an invertible matrix (here (F3̂)); Lemma 0 (here Lemma R2; the extra term
g′(X̂_k) − g′(X̂_0) is bounded); k ≤ T_k ≤ 2k; (O1), (O2). The limit of the heights is Ĵ⁻¹·(Brownian motion with
covariance H) = W by Lemma R2. ∎

A referee should reread the proof of 632 Lemma 2 (c) with these substitutions; nothing else in this file is
"by the same argument".

## 3. The pairing at half time and Jensen's inequality

**Lemma 3.1 (reversal of a set of paths).** For every m and every set Γ of bulk paths (x_0, …, x_m),
P_π((X_0, …, X_m) ∈ Γ) = P̂_π((X̂_m, …, X̂_0) ∈ Γ). In particular, for every bulk state s with floor z and every N,

    P_π(X_0 ∈ cyl(t); heights ≥ 0 during N points; for some k: T_k = N, heights at k equal z, X_k ∈ C(s, z)) = π(C(s, z))·ê_N(s, z).   (3.1)

*Proof.* π(x_0)·ΠP(x_(i−1), x_i) = π(x_m)·ΠP̂(x_i, x_(i−1)). For (3.1): on the left the letters of X_k at depths
z_x + 1, z_y + 1 are the tops −, + of X_0, so z is a floor of X_k; read backward, the path starts in C(s, z), keeps
its heights ≥ −z (the floor letters are not popped), and ends with heights −z at level 0, hence in cyl(t). The
reversed moves read only letters above the floor and the level. ∎

`pairing_check.py` confirms (3.1) on the true chain for every state reachable in N ≤ 7 points (differences below
10⁻¹⁷), with the reversed table of Lemma R0.

**Proposition 3.2 (pairing).** For all n:  R(2n)·π(cyl(t)) ≥ Σ_z ∫ π(ds)·1[z is a floor of s]·ê_n(s, z)·e_n(s, z).

*Proof.* cyl(t) ∩ E_(2n) contains the disjoint events "for some k: T_k = n, X_k = s with heights z; the heights stay
≥ 0 up to k, stay ≥ 0 afterward and return to 0 at level 0 after n more points". Markov property at k, then
(3.1). (The paths in which a stay covers the points n, n + 1 are dropped; `pairing_check.log` shows the size of the
loss for n ≤ 7.) ∎

Fix 0 < a_b < b_b (the box constants; a_b = 1, b_b = 4 below) and box(r) := [a_b·r, b_b·r]² ∩ Z². Let U_n be the following
probability measure on pairs (s, z): z is uniform on box(√n), and s has the law π conditioned on {z is a floor
of s, level(s) ≡ n mod 2} (an event of π-probability 1/8).

**Corollary 3.3 (Jensen).**  R(2n)·π(cyl(t)) ≥ (#box(√n)/8)·exp( −E_(U_n)[−log e_n] − E_(U_n)[−log ê_n] ).

*Proof.* Restrict the sum of Proposition 3.2 to z ∈ box(√n) and to the level of the right parity; it equals
(#box/8)·E_(U_n)[ê_n·e_n], and E[exp(log ê + log e)] ≥ exp(E log ê + E log e). (e_n, ê_n > 0 on the support: (4.1).) ∎

So the lower bound follows from the *typical* descent probability, in logarithmic mean, from a stationary
state with a uniform floor, for the forward and for the reversed chain. No statement about all regular states,
and no statement about the states reached by the growth, is needed.

## 4. The typical descent probability

### 4.0 Three stage bounds and a fallback bound

Box constants a_b = 1, b_b = 4; box(r) = [r, 4r]² ∩ Z², box⁺(r) = [r − 2, 4r + 2]² ∩ Z². A *stage of D points at scale
r from (s, z)* is the run of the chain started at s up to the first move at which the number of points is ≥ D.
"The chain" is P or P̂; every statement of 4.0 holds for both (for P̂ by Lemma R3 in place of 632 Lemma 2).

**(St1) Sharp descent stage** (638 Theorem A, first lines of the proof; uses 632 Lemma 2 (c) and (BM4)). There is
c₄ > 0, depending on the box constants only, such that for every δ ∈ (0, 1/8) there are ε_r, η₀, r₀ with: for r ≥ r₀,
s ∈ A(r; ε_r, η₀, 1), z ∈ box(r) and D with abs(D − (1−δ)r²) ≤ 3: P_s(heights ≥ −z during the stage; z + (heights at its
end) ∈ box(√δ·r)) ≥ q := c₄·δ^(1+μ)/2. (The proof of (BM4) gives the same bound, with the same c₄, for all times
u ∈ [3/4, 1] in place of 1 − δ and target box √δ·[a_b + θ, b_b − θ]²: only u/2 ∈ [3/8, 1/2] is used.)

**(St2) Constant descent stage** (632 Theorem Lo, stage bound; refereed in 639 section 5.6; there the box is
[r, 2r]², and the proof is the same for [r, 4r]²: q₀ > 0 is the infimum of a positive continuous function on a
compact set, the margin is ε′ = 3√δ_*/8). With δ_* := 1/16 there
are absolute constants q_* > 0, ε_*, η_*, r_* with: for r ≥ r_*, s ∈ A*(r) := A(r; ε_*, η_*, 1), z ∈ box(r) and
abs(D − (1−δ_*)r²) ≤ 3: P_s(heights ≥ −z during the stage; floor at its end in box(r/4)) ≥ q_*.

**(St3) Growth stage, upper bound** (632 section 2.1 "one stage, near start" with 639 Proposition B1). For every
ε₁ > 0 there is C₁ such that for every δ ∈ (0, 1/8) there are ε_g, η_g, r_g with: for r ≥ r_g, s ∈ A(r; ε_g, η_g, 1/δ),
z ∈ box⁺(r) and D ≥ (1/δ − 1)r² − 3: P_s(heights ≥ −z during the stage) ≤ q_S := C₁·δ^(μ−ε₁/2).
(Lemma 2 (c) with κ = 1/δ and accuracy δ^(1+μ); Proposition B1 with μ′ = μ − ε₁/2 at a start of norm ≤ 6 + δ.)

**(4.1) Fallback.** There are constants C_f, B such that for every bulk state s, every floor z of s with
z_x, z_y ≤ 4r + 2 and every N ≡ level(s) mod 2 with r² − 1 ≤ N ≤ r² + 1, r ≥ C_f:
−log e_N(s, z) ≤ β(r) and −log ê_N(s, z) ≤ β(r), where β(r) := C_f·r + 2B·log r.

*Proof.* By (F4) the true chain goes from the state "top z letters of s above t" to t in L ≤ 2(z_x + z_y + 1) + 2
points along a fixed path of probability ≥ δ₀^L, and L ≡ level(s) ≡ N mod 2. Then a loop t → t of N − L points:
e_N(s, z) ≥ δ₀^L·R(N − L) and R(2j) ≥ j^(−B) for j ≥ 1 (632 Theorem Lo), R(0) = 1. For ê: by (3.1) and (F4) (from t
every state is reached in L′ ≤ 6(z_x + z_y + 1) + 1 points along a path of probability ≥ δ₀^(L′)),
π(C(s, z))·ê_N(s, z) ≥ π(cyl(t))·R(N − L′)·δ₀^(L′), and π(C(s, z)) ≤ 1. ∎

### 4.1 The theorem

**Theorem 4.1 (typical descent; proved, not reviewed).** For every ε > 0 there is n₀ such that for all n ≥ n₀

    E_(U_n)[−log e_n(s, z)] ≤ (1 + μ + ε)·log n     and     E_(U_n)[−log ê_n(s, z)] ≤ (1 + μ + ε)·log n.

*Proof for e_n (the forward chain descends; the reversed chain controls the density).*

*Parameters.* Absolute: δ_*, q_*, ε_*, η_*, r_*, C_f, B, c₄. Given ε choose θ ∈ (0, 1/4) and ε₁ > 0, then C₁ (St3), then
δ ∈ (0, 1/8), then the constants of (St1) for P and of (St3) for P̂, so that, with

    ϰ := log(2C₁/c₄)/log(1/δ) + ε₁/2 + 2θ·log(1/q_*)/log(16),

    log(2/c₄)/log(1/δ) + 2θ·log(1/q_*)/log(16) + 32·C_f·ϰ/ε_*² ≤ ε/2.

*Scales.* r_k := √n·δ^(k/2) for 0 ≤ k ≤ k_θ, k_θ the largest k with r_k ≥ n^θ; r_k := r_(k_θ)·4^(−(k−k_θ)) for k_θ < k ≤ J, J the
largest k with r_k ≥ r^♯ := (8ϰ/ε_*²)·log n. So n^θ ≤ r_(k_θ) < n^θ·δ^(−1/2), r^♯ ≤ r_J < 4r^♯, k_θ ≤ log n/log(1/δ) and
J − k_θ ≤ (2θ·log n + log(1/δ))/log 16.
*Stages.* The chain starts at (s, z) with n points to go. Stage k (1 ≤ k ≤ J) ends at the first move κ_k at which
the number of points used is ≥ n − ⌊r_k²⌋ (absolute targets; the overshoot is 0 or 1 and is handled as in 639
section 5.1; it is suppressed in the notation). F_k is the event that the floor is not popped during stage k and
that the floor at κ_k lies in box(r_k). Classes: A_k := A(r_k; ε_r, η₀, 1) for k < k_θ and A_k := A*(r_k) for k_θ ≤ k ≤ J;
q_k := q for k ≤ k_θ and q_k := q_* for k > k_θ. By (St1), (St2): for s ∈ A_(k−1) and floor in box(r_(k−1)), P_s(F_k) ≥ q_k.
V_k(s, z) := −log e_(N_k)(s, z), N_k the number of points to go at κ_k.

*Step 1 (Jensen at one stage).* For (s, z) with s ∈ A_(k−1), z ∈ box(r_(k−1)): by the Markov property at κ_k,
e_(N_(k−1))(s, z) ≥ E_s[1_(F_k)·e_(N_k)(X_(κ_k), z′)] = P_s(F_k)·E_s[e_(N_k)(X_(κ_k), z′) | F_k], z′ the floor at κ_k, hence

    V_(k−1)(s, z) ≤ log(1/q_k) + E_s[ V_k(X_(κ_k), z′) | F_k ].

For every (s, z) with z ∈ box(r_k): V_k(s, z) ≤ β(r_k) by (4.1).
*Step 2 (unrolling).* Let P̃_(s,z) be the law of the Markov chain on pairs (state, floor) at the stage ends whose
k-th step is the conditional law of (X_(κ_k), z′) given F_k, and let ϑ be the first k ≥ 0 with X_(κ_k) ∉ A_k, or J if
there is none. Apply Step 1 while k < ϑ and the fallback at ϑ:

    V_0(s, z) ≤ Σ_(k≤J) log(1/q_k) + β(r_J) + Σ_(k<J) β(r_k)·P̃_(s,z)(ϑ = k).

For 1 ≤ k < J, since each conditioning divides by P(F_j) ≥ q_j at a state of A_(j−1):
P̃_(s,z)(ϑ = k) ≤ (q_1 ⋯ q_k)^(−1)·P_s(F_1 ∩ … ∩ F_k ∩ {X_(κ_k) ∉ A_k}).
*Step 3 (the density, by reversal).* U_n has density at most 8/#box(√n) with respect to π ⊗ (counting measure on
box(√n)). **Claim.** For k_θ < k ≤ J and every set B of bulk states,

    Σ_(z∈box(√n)) P_π(F_1(z) ∩ … ∩ F_k(z) ∩ {X_(κ_k) ∈ B}) ≤ #box(r_k)·[ π(B)·q_S^(k_θ) + k_θ·(n+1)·max_(j≤k_θ) π(Â^g_j c) ],

Â^g_j := Â(r_j; ε_g, η_g, 1/δ) (reversed classes). *Proof of the claim.* By Lemma 3.1 the left side is the sum over
z′ ∈ box(r_k) of P̂_π(X̂_0 ∈ B, Γ̂(z′)), where Γ̂(z′) is the reversed image: the reversed chain keeps its heights ≥ −z′,
and at the reversed image of κ_j its floor is in box(r_j), j = k − 1, …, 0. With D_j := ⌊r_j²⌋ − ⌊r_k²⌋, the reversed
image of κ_j has between D_j − 1 and D_j + 1 reversed points. Let ς_j be the first reversed move with at least
D_j − 1 points; the image of κ_j is at most two moves later, so at ς_j the floor is in box⁺(r_j), and between ς_j
and ς_(j−1) there are at least (1/δ − 1)r_j² − 3 points (j ≤ k_θ). Drop all conditions before ς_(k_θ). On the event that
X̂_(ς_j) ∈ Â^g_j for all 1 ≤ j ≤ k_θ, the Markov property of P̂ at ς_(k_θ), …, ς_1 and (St3) for P̂ give the factor π(B)·q_S^(k_θ)
(the first factor is P̂_π(X̂_0 ∈ B) = π(B)). The complement has probability at most Σ_j (n+1)·π(Â^g_j c) by
Lemma R3 (1̂). ∎
*Step 4 (the sums).* (i) k_θ < k < J. By Steps 2, 3 with B = A_k^c, #box(r_k)/#box(√n) ≤ 2r_k²/n = 2δ^(k_θ)·16^(−(k−k_θ)):

    E_(U_n) P̃(ϑ = k) ≤ 16·(δ·q_S/q)^(k_θ)·(1/(16q_*))^(k−k_θ)·π(A*(r_k)^c) + 16·q^(−k_θ)·q_*^(−(k−k_θ))·k_θ·(n+1)·max_j π(Â^g_j c).

δ·q_S/q = (2C₁/c₄)·δ^(−ε₁/2), so the first coefficient is at most C_δ·n^ϰ. By 632 Lemma 1, π(A*(r)^c) ≤ C·r³·exp(−ε_*²r/4), and
r_k ≥ r^♯ gives n^ϰ·exp(−ε_*²r_k/4) ≤ exp(−ε_*²r_k/8). The second term is at most n^C·exp(−ε_g²·n^θ/4) (Lemma R3 (1̂),
r_j ≥ n^θ). So Σ_(k_θ<k<J) β(r_k)·E_(U_n)P̃(ϑ = k) ≤ C_δ·Σ_k r_k⁴·exp(−ε_*²r_k/8) + o(1) = o(1) (the r_k are geometric, r_k ≥ r^♯ → ∞).
(ii) 1 ≤ k ≤ k_θ: P_π(F_1(z) ∩ … ∩ {X_(κ_k) ∉ A_k}) ≤ (n+1)·π(A_k^c) for each z (639 Lemma 1′), so
E_(U_n)P̃(ϑ = k) ≤ 8·q^(−k)·(n+1)·π(A_k^c) ≤ n^C·exp(−c·n^θ), and β(r_k) ≤ C·√n. (iii) k = 0: U_n(s ∉ A_0) ≤ 8π(A_0^c).
(iv) The constants: Σ_(k≤J) log(1/q_k) ≤ k_θ·log(2/(c₄δ^(1+μ))) + (J − k_θ)·log(1/q_*)
≤ (1 + μ + log(2/c₄)/log(1/δ) + 2θ·log(1/q_*)/log 16)·log n + C_δ, and β(r_J) ≤ 4C_f·r^♯ + 2B·log(4r^♯) = (32C_fϰ/ε_*²)·log n + O(log log n).
Together: E_(U_n)[V_0] ≤ (1 + μ + ε/2)·log n + C_δ + O(log log n) ≤ (1 + μ + ε)·log n for n large.

*Proof for ê_n.* Exchange the roles of P and P̂: the reversed chain descends ((St1), (St2), (4.1) for P̂, Lemma R3),
and in Step 3 the reversal of the reversed chain is the forward chain, for which (St3) is 632 section 2.1 and
the exceptional sets are bounded by 632 Lemma 1 and 639 Lemma 1′. ∎

Remarks. (1) *Where the obstruction of 638 section 1 went.* There an irregularity at scale r had to be negligible
against n^(−α), which forces r ≥ C·log n/ε² with the accuracy ε needed for the sharp stage, and the brute-force end
then costs n^(−C′/ε). Here an irregularity at scale r only has to be negligible against n^(−ϰ), the density of the
conditioned law with respect to π, and ϰ → 0 while the accuracy ε_* of the last scales is an absolute constant:
below n^θ the stages are not sharp, and they do not have to be, because there are only θ·log n of them.
(2) The density bound of Step 3 is the place where the reversed chain is needed with its sharp exponent: the
conditioned law at the stage end k is π weighted by a *growth* probability of P̂ from scale r_k to √n, of order
(r_k²/n)^μ, and q^(k_θ) ≈ (r_k²/n)^(1+μ), #box(r_k)/#box(√n) ≈ r_k²/n: the three cancel up to n^ϰ. A forward bound for
the same quantity is a statement of the kind (DR): the reversal replaces it.
(3) Proposition G of 638 (bottom regularity of the growth) is not needed, because the growth is never looked at:
both halves of the excursion are descents, from a stationary middle.

## 5. The theorem

**Theorem 5.1 (proved, not reviewed).** For the meanders with winding function in [−1, 2]:
R(2n) = n^(−α+o(1)) and M_m = Λ^m·m^(−α+o(1)), Λ = (1+2m₀)(1+m₀)² = 10.4245356612…, α = 1 + π/θ = 3.7757713065….

*Proof.* Upper bound: 632 Theorem U (refereed in 639). Lower bound: Corollary 3.3 and Theorem 4.1 give, for every
ε > 0 and n large, R(2n) ≥ (#box(√n)/(8π(cyl t)))·n^(−2(1+μ+ε)) ≥ c·n^(−1−2μ−2ε). M_m = Λ^(m−1)·R(2m − 2) by (F1). ∎

**Dependency list.** Refereed (630, 639): the chain and (F1) to (F4); 625 Theorem 1 (π); 632 Lemmas 0, 1, 2 with 639
Lemma 1′ and sections 3, 5; 632 Theorem U; 632 Theorem Lo; 632 section 2.1 stage bound with 639 Proposition B1;
639 Proposition B2. Not refereed: (BM4) of 638 section 2 (a half-page consequence of B2 and of the positivity and
continuity of the killed heat kernel; used in (St1)); everything in sections 2 to 4 of this file, in particular
Lemma R3 (the transfer of 632 Lemma 2 to P̂; its inputs R0, R1, R2 are exact computations) and Theorem 4.1.
All statements about M_m are conditional, as in 625, on the label machine AE1 and on Propositions 2, 3 of 622 as
refereed in 630.
Not used: (DR*), (DR), (BR), (M), (E); Propositions G, T, D, Lemma 2′, Theorems A, B of 638; Theorem S of 632.

**Recheck of Theorem B of 638 (task 3).** Theorem B derives the same lower bound from (DR*). With Theorem 5.1 it is
not needed. I did not find an error in it; (DR*) stays a conjecture and is now a statement of independent
interest only (regularity of the conditioned chain for every regular start).

**What a referee should attack first.** (1) Lemma R3: reread 632 Lemma 2 (c) and 639 section 3 for P̂. (2) Step 2 and
Step 3 of Theorem 4.1: the bound of P̃(ϑ = k) and the reversed stopping times ς_j. (3) (St1): the uniformity of
(BM4) in the time and its constant c₄ independent of δ (the cancellation in Step 4 needs q ≥ c·δ^(1+μ) with c
independent of δ, and q_S ≤ C·δ^(μ−ε₁/2) with C independent of δ). (4) (4.1) and the parity bookkeeping.
