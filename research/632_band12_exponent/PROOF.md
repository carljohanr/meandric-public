# 632: the band [−1, 2]. Survival exponent, the upper bound M_m ≤ Λ^m·m^(−α+o(1)), and a polynomial lower bound (2026-10-10)

Status words as in `WAVE_620_COMMON.md`. Every theorem below has a complete argument in this file. None of it has
been refereed. What is used from outside is listed in section 0.3; two classical inputs are cited from memory.

## 0. Setting, notation, inputs

### 0.1 Objects (from 625 `PROOF.md` section 0 and 630 `REPAIRS.md`)

m₀ is the positive root of 4x⁴ − 2x − 1; Λ = (1+2m₀)(1+m₀)². M_m is the number of meanders of order m with winding
function in [−1, 2]. The state of the critical chain is (level, σ, τ), σ and τ words over {+, −} written top first.
The *bulk chain* is the move table of 625 section 0 on bottomless stacks; P_s is its law started at the state s; π is
the invariant product law (625 Theorem 1; hand proof in 630 section 4.1); P_π is the stationary bulk chain. Time is
counted in moves (index k) or in points; T_k is the number of points of the first k moves and k(u) is the first k with
T_k ≥ u. z_k = (x_k, y_k) are the heights of σ and τ relative to time 0. δ₀ = 0.088… is the smallest move probability.

t is the start of the true chain; cyl(t) = {level 0, top of σ is −, top of τ is +}; π(cyl(t)) > 0. Q = [0, ∞)² is the
quarter plane. S(u) is the event that x_k ≥ 0 and y_k ≥ 0 for all k ≤ k(u) (the two time-0 tops are not popped during
u points). E_n is the event of 625 (2.1). The facts used:

- (F1) M_m·Λ^(−(m−1)) = P_π(E_(2m−2) | cyl(t)), and on cyl(t) ∩ S(u) the bulk chain is the true chain started at t, up to
  u points (625 (2.1); 630 section 3.3).
- (F2) L¹, L² are martingales of the bulk chain from every state, with increments bounded by B = 3 (630 section 4.2).
- (F3) (625 Lemma 3 and (3.1)) L_k = J·z_k + (1/2)δ^σ·D^σ_k + (1/2)δ^τ·D^τ_k + (a term bounded by 1), det J ≠ 0, with
  D^σ_k = S(σ at time k; b_k) − S(σ at time 0; a_k), a_k = −min_(j≤k) x_j, b_k = x_k + a_k, S(w; i) the sum of the top i
  letters of w. The same for τ.
- (F4) From t every state with heights (x, y) is reached in at most 6(x + y + 1) + 1 points, and t is reached from it in
  at most 2(x + y + 1) + 2 points, all moves having probability ≥ δ₀ (630 section 4.4).

Constants: Σ is the covariance of 625 Theorem 3, θ = arccos((1+2m₀)/(3+4m₀)), **μ = π/(2θ) = 1.3878856532…**,
α = 1 + 2μ = 3.7757713064…. C, c denote positive constants that may change from line to line; they never depend on n.

### 0.2 The quadratic variation in natural time (replaces 625 Proposition 2 (b), with a proof)

**Lemma 0.** There are a symmetric matrix H and constants c_ij ∈ R⁴, γ_ij such that, from every starting state,
A^(ij)_k := Σ_(l<k) E[ΔLⁱ_(l+1)·ΔLʲ_(l+1) | F_l] = H_ij·T_k + c_ij·N(k) + γ_ij·Δ_k + R^(ij)_k, where N(k) is the vector of the four
net letter counts, Δ_k the change of the indicator of level 0, and R^(ij) a martingale with bounded increments.
J⁻¹·H·J⁻ᵀ = Σ.

*Proof.* The conditional law of a move depends on the state a = (level, active top) only: four states. Let
ω = ((1−β)/2, β/2, (1−β)/2, β/2), β = m₀/(1+2m₀), be the π-weights of the four states. The five vectors
(E[ΔN_c | a])_a, c = σ+, σ−, τ+, τ−, and (E[ΔΔ | a])_a are orthogonal to ω (zero stationary fluxes), and they span a
space of dimension 3 (`checks.py`: a 3×3 minor equals 0.119), so they span ω^⊥. ω·(E[ΔT | a])_a = (2+2m₀)/(1+2m₀) ≠ 0.
Hence every function f on the four states is f = h·E[ΔT | ·] + Σ_c c^c·E[ΔN_c | ·] + γ·E[ΔΔ | ·] with
h = ω·f / ω·E[ΔT | ·]. Apply this to f = E[ΔLⁱΔLʲ | ·] and sum along the path; R collects the differences between ΔT,
ΔN, ΔΔ and their conditional expectations. H_ij = E_π[ΔLⁱΔLʲ]/E_π[ΔT], and `checks.py` (40 digits) confirms that
J⁻¹HJ⁻ᵀ has the entries of Σ. ∎

### 0.3 Outside results used

- (O1) Azuma–Hoeffding with Doob: for a martingale M with M₀ = 0 and increments bounded by B,
  P(max_(j≤k) |M_j| ≥ a) ≤ 2·exp(−a²/(2kB²)). Hoeffding's inequality for sums of independent fair signs.
- (O2) The martingale functional central limit theorem (Ethier and Kurtz, Markov Processes, Theorem 7.1.4):
  martingales Mⁿ with jumps tending to 0 uniformly and predictable quadratic variations Aⁿ(u) → H·u in probability
  for each u, with jumps of Aⁿ tending to 0, converge in law to a Brownian motion with covariance H. Cited from
  memory, as in 625.
- (O3) Brownian motion W with covariance Σ in the quarter plane Q (equivalently standard planar Brownian motion in a
  wedge of angle θ); τ_Q is the exit time.
  (BM1) P_ξ(τ_Q > u) ≤ C·(|ξ|/√u)^(2μ) for all ξ ∈ Q, u > 0.
  (BM2) For 0 < a < b and a compact box D inside the open quarter plane there is c > 0 with
  P_ξ(τ_Q > u, W_u ∈ D) ≥ c·|ξ|^(2μ) for all u ∈ [1/2, 1] and all ξ = λζ, ζ ∈ [a, b]², 0 < λ ≤ 1.
  Both follow from the classical series for the heat kernel of a wedge (Carslaw and Jaeger; Spitzer 1958; DeBlassie
  1987; Bañuelos and Smits 1997, Brownian motion in cones): p_u(ξ, y) = (2/(θu))·exp(−(r² + r′²)/(2u))·
  Σ_(j≥1) I_(jπ/θ)(rr′/u)·sin(jπφ/θ)·sin(jπφ′/θ) in polar coordinates of the wedge, whose leading term for small r is of
  order r^(π/θ) = r^(2μ). Cited from memory; not checked at the sources. Gwynne, Mao and Sun (arXiv 1502.00546, Lemma
  3.2) use the same estimates.
- Not used: Denisov–Wachtel, Theorem 2 of 227, any local limit theorem, time reversal, 625 Theorem 3 and its mixing
  form (the uniform Lemma 2 below replaces them).

**Lemma BM3.** sup over ξ ∈ Q of P_ξ(τ_Q > u, |W_u| ≤ R) ≤ C·R^(2+2μ)·u^(−1−μ) for all u > 0, R > 0.

*Proof.* In the linear coordinates in which W is standard, the killed kernel p^Q is symmetric. By Chapman–Kolmogorov,
p^Q_u(ξ, y) = ∫ p^Q_(u/2)(ξ, w)·p^Q_(u/2)(y, w) dw ≤ sup_w p^free_(u/2)(ξ, w)·P_y(τ_Q > u/2) ≤ (C/u)·(|y|/√u)^(2μ) by (BM1). Integrate
over |y| ≤ R. ∎

## 1. Regular states

Fix ε > 0, η > 0, κ > 0 and a scale r ≥ 1. A word w is *ε-regular at scale r* if |S(w; i)| ≤ ε·max(i, r) for all
i ≥ 1. For the bulk chain let G(r, ε, κ) be the event that at every move k ≤ κr² both words σ, τ at time k are ε-regular
at scale r. Define the class

    A(r) = A(r; ε, η, κ) = { s : both words of s are ε-regular at scale r, and P_s(G(r, ε, κ) fails) ≤ η }.

**Lemma 1 (irregular states are exponentially rare).** π(A(r; ε, η, κ)^c) ≤ C_ε·(1 + κr²)·r·exp(−ε²r/4)·(1 + 1/η) for
εr ≥ 4. Consequently, for every event F of the stationary bulk chain and every fixed time k,
P_π(F ∩ {X_k ∉ A(r)}) ≤ π(A(r)^c).

*Proof.* Under π a word is a top letter followed by independent fair coins, so by Hoeffding
P(|S(w; i)| > ε·max(i, r)) ≤ 2·exp(−(ε·max(i, r) − 1)²/(2i)), which is at most 2·exp(−ε²r/4) for i ≤ r and at most
2·exp(−ε²i/4) for i > r. The sum over i is at most C_ε·r·exp(−ε²r/4). By stationarity the same bound holds at each
move k, so P_π(G fails) ≤ 2(1 + κr²)·C_ε·r·exp(−ε²r/4), and π{s : P_s(G fails) > η} ≤ P_π(G fails)/η by Markov's
inequality. The last sentence is stationarity: P_π(X_k ∉ A) = π(A^c). ∎

This is the device that replaces the i.i.d. structure of the hamburger–cheeseburger word: a property of the current
state that fails on a set of π-measure q costs at most q in absolute probability, at every time, whatever the past.

**Lemma 2 (heights from regular states).** Let ε ≤ ε₀ (a constant depending on J only) and s ∈ A(r; ε, η, κ).
(a) On G(r, ε, κ), for all k ≤ κr²: with Z*_k = max_(j≤k) |z_j| and L*_k = max_(j≤k) |L_j|,
    Z*_k ≤ C·L*_k + C·εr + C   and   |z_k − J⁻¹L_k| ≤ C·ε·(r + L*_k) + C.
(b) (concentration) P_s(max over k ≤ k(κr²) of |z_k| ≥ λr) ≤ 4·exp(−cλ²/κ) + η for all λ ≥ 1, r ≥ r₀.
(c) (uniform invariance principle) For every κ and ε′ > 0 there are ε, η₀ > 0 and r₀ such that for all r ≥ r₀ and all
s ∈ A(r; ε, η₀, κ) the law under P_s of the path u ↦ r⁻¹·z_(k(r²u)), 0 ≤ u ≤ κ (linearly interpolated between integer
values of r²u), is within Prokhorov distance ε′, for the uniform metric on C[0, κ], of the law of W started at 0.

*Proof.* (a) By (F3), |D^σ_k| ≤ ε·max(b_k, r) + ε·max(a_k, r) ≤ ε·(2r + 3Z*_k) on G (the word at time 0 is regular because
s ∈ A; the word at time k because G holds), and the same for τ. So |z_k − J⁻¹L_k| ≤ C·ε·(r + Z*_k) + C, hence
Z*_k ≤ ‖J⁻¹‖·L*_k + C·ε·(r + Z*_k) + C, and for C·ε ≤ 1/2 the first claim follows; insert it into the previous line.
(b) By (a), on G the event Z* ≥ λr forces L* ≥ cλr for λ ≥ 1, ε ≤ ε₀, r ≥ r₀. The number of moves is at most κr². Use
(O1) for L¹ and L², and P_s(G^c) ≤ η.
(c) Suppose not. Then there are ε′ > 0, r_n → ∞, ε_n → 0, η_n → 0 and s_n ∈ A(r_n; ε_n, η_n, κ) with distance > ε′. Under
P_(s_n) put Mⁿ(u) = r_n⁻¹·L_(k(⌊r_n²u⌋)). k(·) is a bounded stopping time, so Mⁿ is a martingale with jumps ≤ B/r_n, and
Aⁿ(u) = r_n⁻²·A_(k(⌊r_n²u⌋)) (Lemma 0) is its predictable quadratic variation, with jumps ≤ B²/r_n². By Lemma 0,
Aⁿ(u) − H·u = r_n⁻²·(c·N + γ·Δ + R) + O(r_n⁻²). Here |N(k)| ≤ |z_k| + |D^σ_k| + |D^τ_k| ≤ C·(r_n + Z*) on G, Z*/r_n is tight by
(b), and max_k |R_k| = O_P(r_n) by Doob. So Aⁿ(u) → H·u in probability, and by (O2) Mⁿ converges to a Brownian motion
with covariance H. By (a), sup_u |r_n⁻¹·z_(k(r_n²u)) − J⁻¹Mⁿ(u)| ≤ C·ε_n·(1 + L*/r_n) + C/r_n on G, which tends to 0 in
probability. So the height paths converge in law to J⁻¹·(Brownian motion with covariance H) = W; the interpolation
changes the path by at most C/r_n. This contradicts distance > ε′. ∎

Use of (c): for a closed set F of paths, P_s(path ∈ F) ≤ P(W ∈ F^(ε′)) + ε′ and P_s(path ∈ F^(ε′)) ≥ P(W ∈ F) − ε′, where
F^(ε′) is the ε′-neighborhood. Survival events are of this kind: a discrete path that stays in Q − z has its
interpolation in Q − z.

Remark (what is not proved). Lemma 2 is about states of A(r). That a state with regular letters automatically
satisfies the dynamic condition P_s(G fails) ≤ η is not proved; Lemma 1 only says that the exceptions are rare under π.
This is the source of the gap in section 6.

## 2. The survival exponent (weak form)

**Theorem S.** P_π(S(n) | cyl(t)) = n^(−μ+o(1)) and P_π(S(n)) ≤ n^(−μ+o(1)), μ = π/(2θ). In words: the probability that the
true chain started at t is not killed during n points is n^(−π/(2θ)+o(1)). Moreover, for 0 < a < b suitably chosen, the
lower bound holds with the additional requirement that the heights after n points lie in [a√n, b√n]².

Fix ε₁ ∈ (0, 1/4). δ ∈ (0, 1/8) is a small parameter (chosen at the end, depending on ε₁). J = ⌊log n/log(1/δ)⌋.
Scales: r_j = √n·δ^((J−j)/2) for 0 ≤ j ≤ J, so r_J = √n and 1 ≤ r_0 < δ^(−1/2); times T_j = r_j² (integer parts are omitted; a stage ends at the first move at which the number of points reaches
its target, which shifts every time by at most 1 and changes nothing below).

### 2.1 Upper bound

Let ρ = δ^(μ−ε₁), K = K(δ) with c·K²/8 = 3·log(1/ρ) (c the constant of Lemma 2 (b) for κ = 1/δ; since κδ < 1 the bound of
Lemma 2 (b) for displacements measured in units of r/√δ is 4·exp(−cλ²/4) + η), ε′ = δ^(1+μ), and ε, η₀, r₀ from Lemma 2 (c)
for κ = 1/δ. Put η = n^(−10) and A_j = A(r_j; ε, η, 1/δ). Fix ε₂ ∈ (0, 1) and j₁ = ⌈ε₂J⌉, so r_(j₁) ≥ n^(ε₂/2)/C. Let
S′_j = {z_k ∈ Q for k(T_(j₁)) ≤ k ≤ k(T_j)}, Reg_j = {X_0 ∈ A(r_(j₁); ε, η, 1), X at time T_i in A_i for j₁ ≤ i ≤ j},
p_j = P_π(S′_j ∩ Reg_j), f_j(λ) = P_π(S′_j ∩ Reg_j ∩ {|z at time T_j| > λr_j}).

*One stage, near start.* For s ∈ A_j and |z| ≤ Kr_j: by Lemma 2 (c) and (BM1),
P_s(z + z_k ∈ Q during (1/δ − 1)·r_j² points) ≤ P_(ζ+ε′(1,1))(τ_Q > 1/δ − 1) + ε′ ≤ C·((K + 1)²·δ)^μ + δ^(1+μ) =: q_S (ζ = z/r_j).
*Recursion.* By the Markov property at time T_j, p_(j+1) ≤ q_S·p_j + f_j(K). If |z at T_j| > λr_j, then either the
displacement during stage j exceeds λr_j/2 = (λ/(2√δ))·r_(j−1), or |z at T_(j−1)| > (λ/(2√δ))·r_(j−1). By Lemma 2 (b),
f_j(λ) ≤ p_(j−1)·(4·exp(−cλ²/4) + η) + f_(j−1)(λ/(2√δ)), and f_(j₁)(λ) ≤ 4·exp(−cλ²) + η (stage from time 0). Since
(4δ)⁻¹ ≥ 2, iteration gives
f_j(K) ≤ Σ_(i=1..j−j₁) p_(j−i)·4·exp(−c·K²·2^(i−1)/4) + 4·exp(−c·K²·2^(j−j₁)/4) + (J + 1)·η.
*Induction.* Claim p_j ≤ 2ρ^(j−j₁) for j₁ ≤ j ≤ J. With the choice of K, exp(−c·K²·2^i/8) ≤ ρ^(3·2^i) ≤ ρ^(2+i+1)·ρ^i (because
3·2^i ≥ 2i + 3), so Σ_i ρ^(−i)·4·exp(−cK²2^(i−1)/4) ≤ 8ρ² and f_j(K) ≤ ρ^(j−j₁)·(16ρ² + 4ρ²) + (J + 1)·n^(−10). Hence
p_(j+1) ≤ 2ρ^(j−j₁)·(q_S + 10ρ²) + (J + 1)·n^(−10). For δ small, q_S ≤ C·log^μ(1/δ)·δ^μ ≤ ρ/4 and 10ρ² ≤ ρ/4, and
(J + 1)·n^(−10) ≤ ρ^(J) ≤ ρ^(j+1−j₁) for n large (ρ^J ≥ n^(−μ)). The claim follows.
*Conclusion.* P_π(S(n)) ≤ p_J + P_π(Reg_J fails) ≤ 2ρ^(J−j₁) + (J + 2)·max_i π(A_i^c). By Lemma 1 the last term is at most
C·n^(12)·exp(−ε²·n^(ε₂/2)/C), and ρ^(J−j₁) = (δ^J)^((μ−ε₁)(1−ε₂)+o(1)) with δ^J = r_0²/n ≤ 1/(δn). So
P_π(S(n)) ≤ n^(−(μ−ε₁)(1−ε₂)+o(1)). ε₁, ε₂ are arbitrary. Finally P_π(S(n) | cyl(t)) ≤ P_π(S(n))/π(cyl(t)). ∎

### 2.2 Lower bound

Fix 0 < a < b (for instance a = 1, b = 2), box_j = [a·r_j, b·r_j]², and let
g_j = P_π(cyl(t) ∩ S(T_j) ∩ {z at time T_j ∈ box_j}). Let ε′ = δ^(1+μ), ε, η₀, r₀ from Lemma 2 (c) with κ = 1/δ, and
A_j = A(r_j; ε, η₀, 1/δ) (here η₀ is a constant). c₁ = ε²/4.

*One stage.* For s ∈ A_j and z ∈ box_j, by Lemma 2 (c), scaling by √δ, and (BM2) with D = [1.1a, 0.9b]²:
P_s(z + z_k ∈ Q during the stage, end in box_(j+1)) ≥ P_(√δ(ζ−ε′(1,1)))(τ_Q > 1 − δ, W_(1−δ) ∈ D) − ε′ ≥ c·δ^μ − δ^(1+μ) =: q_L ≥ c·δ^μ/2.
*Recursion.* By the Markov property and Lemma 1, g_(j+1) ≥ q_L·(g_j − π(A_j^c)) ≥ q_L·(g_j − C_δ·r_j³·exp(−c₁r_j)).
*Start.* Claim: there is C₄ such that for every ε_fl ∈ (0, a/2) and all r ≥ r(ε_fl),
g(r) := P_π(cyl(t) ∩ S(r²) ∩ {z at time r² ∈ [ar, br]²}) ≥ exp(−C₄·ε_fl·r).
Proof: let d = ⌈ε_fl·r⌉. For words u, v of length d let s(u, v) be the state of the true chain with σ = u and τ = v on
top of the letters of t. By (F4) there is a fixed path γ(u, v) from t to s(u, v) with at most 12d + 7 points and
probability ≥ δ₀^(12d+7); no γ(u, v) is a prefix of another (the building procedure of 625 builds σ completely while
τ = (+), so the state s(u, v), in which τ has d + 1 letters, occurs only at the end of γ(u, v)). Let F* be the event,
for the bulk chain, that the heights stay ≥ −(d − 1) during r² points and lie in [ar + εr, br − εr]² − (d, d) during
the last 12d + 7 points. F* depends on the level, the uniform variables and the top d letters of each stack only.
The paths "γ(u, v), then F* for the remaining points" are disjoint for different (u, v) and are contained in the
event of g(r). So g(r) ≥ π(cyl(t))·δ₀^(12d+7)·Σ_(u,v) P_π(F* | level 0, top words u, v), and since π(level 0, u, v) ≤ 2^(−2d),
the sum is at least 2^(2d)·P_π(F* ∩ level 0). By Lemma 2 (c) (κ = 1) and Lemma 1, P_π(F* ∩ level 0) ≥ c(ε_fl) > 0 for
r large (a Brownian path stays above a negative level and ends, and stays for a short final interval, in an open
set, with positive probability). This gives the claim with C₄ = 13·log(1/δ₀).
*Induction.* Choose ε_fl = c₁/(4C₄) and j₀ with r_(j₀) large. Then g_(j₀) ≥ exp(−c₁r_(j₀)/4), and by induction
g_j ≥ g_(j₀)·(q_L/2)^(j−j₀) for j ≥ j₀: the induction step needs C_δ·r_j³·exp(−c₁r_j) ≤ g_j/2, which holds because r_j grows
geometrically in j while (q_L/2)^(j−j₀) decays only geometrically (take r_(j₀) ≥ r(δ), a constant depending on δ; the
Lemma 2 (c) threshold r₀ is included in r(δ)). Then j₀ is bounded in n and
g_J ≥ c_δ·(c·δ^μ/4)^(J) = n^(−μ−o_δ(1)+o(1)), where o_δ(1) = log(4/c)/log(1/δ) → 0 as δ → 0. ∎

## 3. The upper bound for the return probability

**Theorem U.** P_π(E_(2n) | cyl(t)) ≤ n^(−1−2μ+o(1)). Hence **M_m ≤ Λ^m·m^(−α+o(1))**, α = 1 + π/θ = 3.7757713064….

The proof follows the heights only. The letters enter only through the classes A(r). No local limit theorem and
no gluing of words is needed: the factor n^(−1) comes from localizing the heights once per scale during the descent.

*Setting.* Fix ε₁ > 0, δ small, ρ = δ^(1+μ−ε₁), K with c·K²/8 = 3·log(1/ρ) (c from Lemma 2 (b) with κ = 1), ε′ = δ^(2+μ),
ε, η₀, r₀ from Lemma 2 (c) with κ = 1. Descent scales r_k = √n·δ^(k/2), k = 0, …, J, where J is the largest k with
r_k ≥ r_* := C_*·log n; C_* is chosen so that, with η = n^(−10) and A_k = A(r_k; ε, η, 1), Lemma 1 gives π(A_k^c) ≤ n^(−12) for
all k ≤ J (possible since π(A_k^c) ≤ C·n^(11)·log³n·exp(−ε²C_*·log n/4)). Stage k (1 ≤ k ≤ J) is the time interval in which the
number of points still to come decreases from r_(k−1)² to r_k²; it has (1 − δ)·r_(k−1)² points.

For a bulk state s ∈ A_k and a floor z ∈ N² (the heights of s above the letters of t) let ũ_k(s, z) be the probability,
for the bulk chain started at s, that the heights stay ≥ −z during r_k² points, equal −z at the end, and that the
states at the ends of the stages k + 1, …, J lie in A_(k+1), …, A_J. Let W̃_k = sup over s ∈ A_k and z ∈ N² of ũ_k(s, z), and
Y_k(λ) the same supremum over |z| > λr_k.

**Claim.** W̃_k ≤ w_k := ρ^(J−k) for 0 ≤ k ≤ J, for n large.

*Proof.* k = J is trivial. Far floors: a return from a floor with |z| > λr_k needs, during stage k + 1 (at most r_k²
points), a displacement ≥ λr_k/2, or else the floor at the end of the stage is still farther than
λr_k/2 = (λ/(2√δ))·r_(k+1). By Lemma 2 (b) and the Markov property, Y_k(λ) ≤ (4·exp(−cλ²/4) + η)·W̃_(k+1) + Y_(k+1)(λ/(2√δ)), and
Y_J(λ) ≤ 4·exp(−cλ²) + η. With the induction hypothesis W̃_j ≤ w_j for j > k, (4δ)⁻¹ ≥ 2 and the choice of K, exactly as in
2.1: Y_k(K) ≤ w_k·Σ_(l≥0) 4·ρ^(−(l+1))·exp(−cK²2^l/4) + (J + 1)·η ≤ 8ρ²·w_k + (J + 1)·η.
Now let s ∈ A_(k−1), z ∈ N². By the Markov property at the end of stage k, with z′ the floor there,
ũ_(k−1)(s, z) ≤ P_s(heights ≥ −z during the stage, |z′| ≤ Kr_k)·W̃_k + Y_k(K).
By Lemma 2 (c) (the event is closed), with ζ = z/r_(k−1), then Lemma BM3 with u = 1 − δ ≥ 1/2 and R = K√δ + 3ε′ ≤ 2K√δ:
P_s(heights ≥ −z, |z′| ≤ Kr_k) ≤ P_(ζ+ε′(1,1))(τ_Q > 1 − δ, |W_(1−δ)| ≤ K√δ + 3ε′) + ε′ ≤ C·(2K)^(2+2μ)·δ^(1+μ) + δ^(2+μ).
So ũ_(k−1)(s, z) ≤ w_k·δ^(1+μ)·(C·log^(1+μ)(1/δ) + δ + 8δ^(1+μ−2ε₁)) + (J + 1)·n^(−10) ≤ w_k·ρ/2 + w_(k−1)/2 = w_(k−1) for δ ≤ δ(ε₁) and n
large (w_(k−1) ≥ ρ^J ≥ n^(−1−μ)). ∎

*Proof of Theorem U.* On cyl(t) ∩ E_(2n): S(n) holds; then, from the state X at time n with floor z (its heights),
the heights stay ≥ −z and return to −z after n more points. Let Bad be the event that the state at time n or at the
end of some stage is outside its class A_k. By Lemma 1, P_π(Bad) ≤ (J + 1)·n^(−12). By the Markov property at time n,
P_π(cyl(t) ∩ E_(2n)) ≤ P_π(S(n))·W̃_0 + P_π(Bad) ≤ n^(−μ+o(1))·ρ^J + n^(−11),
and ρ^J = (δ^J)^(1+μ−ε₁) = (r_J²/n)^(1+μ−ε₁) ≤ (C·log²n/n)^(1+μ−ε₁) (r_J ≤ r_*/√δ). So P_π(E_(2n) | cyl(t)) ≤ n^(−1−2μ+ε₁+o(1)); ε₁ is
arbitrary. With (F1): M_m = Λ^(m−1)·P_π(E_(2m−2) | cyl(t)) ≤ Λ^m·m^(−α+o(1)). ∎

## 4. A polynomial lower bound for the return probability

**Theorem Lo.** There is a constant B such that P_π(E_(2n) | cyl(t)) ≥ n^(−B) for all n ≥ 1. Hence
**M_m ≥ Λ^m·m^(−B)** for all m ≥ 2. B is not explicit.

*Proof.* Fix δ = 1/16, a = 1, b = 2, and the margin ε′ = (b − a)·√δ/8. Let
q₀ = inf over ζ ∈ [a, b]² of P_ζ(W stays in Q + ε′(1,1) during [0, 1 − δ], W_(1−δ) ∈ [a√δ + ε′, b√δ − ε′]²) > 0
(positive and continuous in ζ; elementary). Let ε, η₀, r₀ be as in Lemma 2 (c) for κ = 1 and the accuracy
ε″ = min(ε′, q₀/4). Descent scales r_k = √n·δ^(k/2), k ≤ J, J the largest k with r_k ≥ C_*·log n,
A_k = A(r_k; ε, η₀, 1). Stage event F_k: during stage k the heights stay ≥ 0 (above t) and the floor at its end lies in
box_k = [ar_k, br_k]². By Lemma 2 (c), for s ∈ A_(k−1) and floor in box_(k−1): P_s(F_k) ≥ q₀ − ε″ ≥ q := q₀/2 (the
ε″-neighborhood of the event of q₀ is inside the event F_k).
*Chain.* For s ∈ A_0 with floor in box_0, induction over the stages and the Markov property give
P_s(F_1 ∩ … ∩ F_J) ≥ q^J − Σ_(k<J) P_s(state at the end of stage k ∉ A_k).
*Last scale.* From a state in A_J with floor in box_J: with probability ≥ q′ > 0 (Lemma 2 (c) again: the Brownian
motion from ζ ∈ [a, b]² stays in Q + ε′(1,1) and inside the ball of radius 3b during one unit of time with probability
bounded below) the chain stays above t with both heights at most 4b·r_J during r_J² points, in particular during the
first r_J² − C″·r_J points; then, by (F4), t is reached in at most 16b·r_J + 6
points with probability ≥ δ₀^(16b·r_J+6); the remaining number of points is even (every path from t to t has an even
number of points) and at most C″·r_J for C″ = 16b + 8, and is filled with pairs of moves "go, nothing changes" at t.
This costs at least exp(−C₅·r_J) ≥ n^(−C₅C_*/√δ).
*Growth.* By Theorem S (lower bound with the box; there the parameter δ is taken so small that the exponent is at
most μ + 1), g := P_π(cyl(t) ∩ S(n) ∩ {z at time n ∈ [a√n, b√n]²}) ≥ n^(−μ−1) for n large.
*Assembly.* By the Markov property at time n and Lemma 1 (stationarity at time n and at the end of each stage),
P_π(cyl(t) ∩ E_(2n)) ≥ [ q^J·(g − π(A_0^c)) − Σ_(k≤J) π(A_k^c) ]·q′·n^(−C₅C_*/√δ).
q^J ≥ n^(−b₁) with b₁ = log(1/q)/log(1/δ). By Lemma 1, π(A_k^c) ≤ C·log³n·n·exp(−ε²·C_*·log n/4) for k ≤ J. Choose C_* with
ε²C_*/4 ≥ b₁ + μ + 4. Then the bracket is at least n^(−b₁−μ−1)/2, and the theorem holds with any
B > b₁ + μ + 1 + C₅C_*/√δ for n large; enlarge B for the finitely many remaining n (the probability is positive). ∎

Why this stops at a polynomial with a poor exponent: below the scale C_*·log n the rare-state device of Lemma 1 is
no longer negligible against polynomial probabilities, and the last C_*·log n letters are removed by brute force.

## 5. Summary of what is proved

For all m ≥ 2:   Λ^m·m^(−B) ≤ M_m ≤ Λ^m·m^(−α+o(1)),   α = 1 + π/arccos((1+2m₀)/(3+4m₀)) = 3.7757713064…,
and the true chain survives n points with probability n^(−π/(2θ)+o(1)). In particular lim sup log(M_m·Λ^(−m))/log m = −α is
NOT proved; proved is lim sup ≤ −α and lim inf ≥ −B.

Dependencies: 625 Theorem 1, identity (2.1), Lemma 3 and (3.1); 630 `REPAIRS.md` sections 3 and 4 (hand proofs, (F4));
Lemma 0 here (one rank computed by `checks.py`); (O1), (O2), (O3).

## 6. What is missing for M_m = Λ^m·m^(−α+o(1)): one lemma

**Proposition R (reduction; proved).** Suppose the following holds.
(M) For every small δ there are r₀ and classes Â_r of bulk states, r ≥ r₀, with π(Â_r^c) ≤ exp(−r^c) for some c > 0, such
that for all r ≥ r₀/√δ, all s ∈ Â_r and all floors z ∈ [ar, br]²:
P_s(heights ≥ −z during (1 − δ)r² points; floor at the end in [a√δ·r, b√δ·r]²; state at the end in Â_(√δ·r)) ≥ δ^(1+μ+o_δ(1)).
Then P_π(E_(2n) | cyl(t)) ≥ n^(−α−o(1)), and with Theorem U, M_m = Λ^m·m^(−α+o(1)).

*Proof.* As in section 4, with the stage bound of (M) in place of q, down to the scale r₀ (a constant), where the
remaining removal costs a constant. The state at time n is in Â_(√n) up to an event of probability exp(−n^(c/2)) (Lemma 1's
argument), and P(growth) ≥ n^(−μ−o(1)) by Theorem S. The product over J = log n/log(1/δ) + O(1) stages is
δ^((1+μ+o_δ(1))·J) = n^(−1−μ−o_δ(1)). ∎

The difference between (M) and what section 4 uses: in (M) the end state must be regular for the next scale with a
conditional probability close to 1 given the start state, at every scale down to a constant. Lemma 1 gives this
only in absolute probability, which is enough above the scale log n and not below. (M) would follow from two
statements, neither proved here:
- (DR, dynamic regularity) from a state whose letters are regular in all windows of length ≥ ℓ, the letters pushed
  during the next r² moves and still present are regular in all windows of length ≥ ℓ, up to conditional
  probability exp(−cℓ). Under π this is Lemma 1; the open point is a deterministic regular start (a mixing statement
  for the pair (level, active top) seen from a fixed state).
- (BR, bottom regularity) under P_π(· | cyl(t) ∩ S(n)), the lowest i letters of the two words at time n are regular
  seen from the bottom, up to conditional probability → 0 as i → ∞ uniformly in n.

(M) was not tested numerically: it is an asymptotic statement about conditional laws at scales below log n, and no
exact computation at small size discriminates it. Status of (M), (DR), (BR): **conjecture**.

## 7. Comparison with Gwynne–Mao–Sun (arXiv 1502.00546) and Gwynne–Sun (arXiv 1505.03375)

I read 1505.03375 (sections 1, 5 in full; 2 to 4 for structure). Its Theorem 1.9 is P(X(1, 2n) = ∅) = n^(−1−2μ+o(1)).

| Their ingredient | Role | Here |
|---|---|---|
| i.i.d. symbols; the word after a stopping time is independent of the past | Markov property at every scale, for free | The chain is Markov in the full state (two words). The past enters through the letters. Replaced by the classes A(r) and Lemma 1 (stationarity makes bad states cost their π-measure) |
| Sheffield's Theorem 2.5 (scaling limit of the counts) | Brownian approximation at every scale | Lemma 2 (c): uniform invariance principle from every state of A(r), via the two exact martingales and Lemma 0 |
| Few flexible orders (GMS Lemma 3.7, Corollary 6.2; GS Lemma 5.2) | the mismatch between the counts and the walk is of lower order | Lemma 2 (a): the mismatch is the letter imbalance, at most ε·(scale) on regular states |
| GMS Lemma 3.2, Shimura: cone estimates for Brownian motion | exponents | (BM1), (BM2), Lemma BM3 |
| Lemmas 1.12, 1.13 (regular variation of the no-order time) | first half, exponent μ | Theorem S, in the weaker form n^(−μ+o(1)); regular variation is not proved here |
| Section 5.1 (lower bound): stages [−δ^(k−1)n, −δ^k·n], Lemma 5.1 gives δ^(μ+1) per stage, last finitely many symbols at constant cost | lower bound n^(−1−2μ−o(1)) | Same scheme in sections 4 and 6. It closes only down to the scale C·log n (Theorem Lo). Their step "the same holds if we condition on the past word" has no analog below that scale: this is (M) |
| Sections 2 to 4 (local estimates through the i.i.d. pairs (J^H_m, L^H_m) and Doney's bivariate local limit theorem) and section 5.3 (upper bound through stopping times at scales n^(ξ^k)) | upper bound | Not needed and not available: the times at which one stack reaches a new minimum do not give i.i.d. pieces (the other stack carries the past; 625 AL4). The upper bound here uses geometric stages and the sup recursion of section 3 instead; the factor n^(−1) is the product of the per-stage localization costs of Lemma BM3 |

A remark on the analogy (`checks.py` (a)): the critical chain is driven by i.i.d. symbols of two classes, A with
probability 1/(1+m₀) and B with probability m₀/(1+m₀). A symbol B is a go without pop. A symbol A is a stay if the
active top is + and a go with pop if the active top is −. So A is a "flexible" symbol whose effect depends on the top
letter, as Sheffield's flexible order. The difference from Sheffield's model is that here the top letter decides the
level and the number of points, not only which stack loses a letter, so the reduced effect of a block of symbols
is not a word but a decision tree, and blocks read between record times are not independent.
