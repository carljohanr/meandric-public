# 638: band [−1, 2], the matching lower bound. What is proved, and the one statement that remains (2026-10-10)

Status words as in `WAVE_620_COMMON.md`. **The lower bound M_m ≥ Λ^m·m^(−α−o(1)) is NOT proved here.** Proved here: the
bootstrap of the brief does not close and why (section 1); the descent with the sharp exponent down to the scale
M·log n (Theorem A, section 2); two unconditional statements that remove half of what 632 left open (Proposition G:
bottom regularity (BR) can be built into the growth; Proposition T: the end game costs exp(−o(h)) from all states
outside a set of exponentially small stationary measure; section 3); a comparison inequality between the chain from
an arbitrary state and the stationary chain (Proposition D, section 4); a reduction of the theorem to one
statement about pushed letters, (DR*) (Theorem B, section 4); exact flux identities (Lemma F, section 5).
Everything rests on 632 `PROOF.md` (Lemmas 0, 1, 2, Theorems S, U, Lo; **not refereed**), on 625 Theorem 1 and
Corollary 1, and on 630 `REPAIRS.md` sections 3, 4.

## 0. Notation (all from 632 `PROOF.md` section 0)

m₀, Λ, the bulk chain, π, P_s, P_π, moves and points, T_k, k(u), heights z_k = (x_k, y_k), cyl(t), S(u), E_n, (F1) to
(F4), μ = π/(2θ), α = 1 + 2μ, δ₀ (smallest move probability), the classes A(r; ε, η, κ), Lemmas 0, 1, 2, (O1) to (O3),
(BM1), (BM2), Lemma BM3. W is the Brownian motion with covariance Σ, Q the quarter plane, τ_Q its exit time.
R(2n) := P_π(E_(2n) | cyl(t)); by (F1), M_m = Λ^(m−1)·R(2m − 2).
For a bulk state s and a *floor* z ∈ N² (such that the letters of s at depths z_x + 1 in σ and z_y + 1 in τ are −, +
and the bulk chain restricted to heights ≥ −z is the true chain above t), and N ≡ level(s) mod 2, put

    e_N(s, z) := P_s(heights ≥ −z during N points; level 0 and heights = −z after N points)
               = P(the true chain started at the state "top z letters of s above t" is at t after N points).

Parity: a go has one point and changes the level, a stay has two points, so a path from level l to t has N ≡ l mod 2.
A stage "of u points" ends at the first move at which the number of points is ≥ its target (shift at most 1, as in
632). On the event cyl(t) ∩ (alive) the number of points so far is ≡ level mod 2, so the number of points left
before 2n always has the right parity.

*A remark on Lemma 1 of 632 at stage ends.* Stage ends are stopping times k(T), not fixed move indices. The bound
used below is P_π(F ∩ {X_(k(T)) ∉ A}) ≤ Σ_(k≤T) P_π(X_k ∉ A) ≤ T·π(A^c) (stationarity at each fixed move index). The
extra factor T ≤ 2n is harmless wherever 632 uses Lemma 1 (the bounds there are n^(−12) or stretched exponentials).

## 1. The bootstrap of the brief (idea 3) does not close; the exact obstruction

*The proposal.* Run the stages of 632 down to the scale h = C·log n; replace the brute-force removal of the last
letters (cost exp(−c·h) = n^(−c·C)) by a polynomial bound h^(−C′), obtained by applying Theorem Lo at the scale h.

*Why it fails as posed.* Theorem Lo bounds R(2N) from below: it is a statement about paths that start at t. At the
end of the descent the chain is at a state s with heights of order h, and what is needed is a lower bound for
e_N(s, z), N ≈ h², *for the states s that the conditioned chain actually reaches*. The only information about s is
s ∈ A(h; ε, η, 1): its letter sums are at most ε·max(i, h) at depth i, and the same along the next h² moves with
probability 1 − η. This says nothing about windows shorter than εh. To continue the descent below h one needs the
state at the end of each further stage to be regular at the smaller scale r < h, and Lemma 1 gives this only in
absolute probability exp(−ε²r/4), which is not small compared with the probability n^(−1−2μ) of the path so far. So
the hypotheses of the stage argument inside Theorem Lo are not available below h: the composition does not close.

*The obstruction is quantitative and does not depend on the scheme* (a remark, not a theorem). Any argument that
controls letters through the stationary law pays as follows. An accuracy ε for letter sums at scale r holds up to
stationary probability exp(−ε²r/4) (Hoeffding; this rate is sharp). To be negligible against n^(−α) one needs
r ≥ C·log n/ε². The heights are then known up to ε·r ≥ C·log n/ε letters, and removing that many letters by fixed
moves costs δ₀^(C·log n/ε) = n^(−C′/ε). Every variant I tried reduces to this inequality: thin top layers removed by
brute force, cuts of the growth at a frozen state, classes defined by iterated conditional probabilities (the
recursion of Proposition T below loses its rate at the same place), stages run relative to a frozen floor. In
short: *at the probability level n^(−C) the top C·log n letters of a stack can be arbitrary, and below the scale
log n only statements that hold conditionally on the current state can be used.* The polynomial bound of the
brief, e_N(s, z) ≥ h^(−C′), is needed for states s whose top letters are not controlled; that is statement (E) of
section 2, and it is open.

*What does close* is the other half of the proposal: above the scale M·log n the descent can be run with the sharp
cost per stage. That is Theorem A.

## 2. The descent with the sharp exponent down to the scale M·log n

**(BM4).** For 0 < a < b there is c₄ > 0 such that for all δ ∈ (0, 1/8), all ζ ∈ [a/2, b]² and every θ ≤ a√δ/8:
P_ζ(W stays in Q + θ(1,1) during [0, 1 − δ], W_(1−δ) ∈ [a√δ + θ, b√δ − θ]²) ≥ c₄·δ^(1+μ).

*Proof.* Translate by θ(1,1): the start becomes ζ′ ∈ [a/4, b]², the cone becomes Q and the target contains
√δ·D₁, D₁ = [5a/4, 3b/4]² (take b ≥ 4a, which is allowed). Let p^Q be the killed transition density; it is symmetric
in the coordinates where W is standard, and u = 1 − δ ∈ [7/8, 1]. For y ∈ √δ·D₁ and a compact box D′ inside the
open quarter plane, by Chapman–Kolmogorov and symmetry,
p^Q_u(ζ′, y) = p^Q_u(y, ζ′) ≥ ∫_(D′) p^Q_(u/2)(y, w)·p^Q_(u/2)(w, ζ′) dw ≥ c′·P_y(τ_Q > u/2, W_(u/2) ∈ D′) ≥ c′·c·abs(y)^(2μ),
where c′ > 0 is the infimum of the (positive, continuous) killed density p^Q_v(w, ζ′) over w ∈ D′, ζ′ ∈ [a/4, b]²,
v ∈ [7/16, 1/2], and the last step is (BM2) with time rescaled. abs(y) ≥ a√δ and the area of √δ·D₁ is of order δ. ∎

(BM4) rests on (BM2), which 632 cites from memory, and on positivity and continuity of the killed heat kernel in
the open wedge (classical).

**Theorem A (sharp descent to the logarithmic scale).** For every ε > 0 there are δ ∈ (0, 1/8), regularity parameters
(ε_r, η₀), a constant M and n₀ such that for n ≥ n₀ the following holds. Let r_k = √n·δ^(k/2), let J be the largest k with
r_k ≥ M·log n, box_k = [a·r_k, b·r_k]², A_k = A(r_k; ε_r, η₀, 1), and let F_k be the event that during stage k (the
(1 − δ)·r_(k−1)² points in which the number of points still to come before 2n decreases from r_(k−1)² to r_k²) the heights
stay ≥ 0 and are in box_k at its end. Then

    P_π( cyl(t) ∩ S(n) ∩ {z at time n ∈ box_0} ∩ F_1 ∩ … ∩ F_J ∩ {state at the end of stage J ∈ A_J} ) ≥ n^(−1−2μ−ε).

*Proof.* The proof of Theorem Lo of 632 with one change: the stage probability. Fix δ so small that
κ_δ := log(4/c₄)/log(1/δ) ≤ ε/(4(1+μ)) and that Theorem S (lower bound with the box) holds with exponent μ + ε/4. Let
θ = a√δ/8, ε″ = min(θ, c₄δ^(1+μ)/2), and (ε_r, η₀, r₀) from Lemma 2 (c) of 632 for κ = 1 and accuracy ε″. For s ∈ A_(k−1) and
a floor in box_(k−1), Lemma 2 (c) and (BM4) give P_s(F_k) ≥ c₄δ^(1+μ) − ε″ ≥ q := c₄δ^(1+μ)/2 (the ε″-neighborhood of the
Brownian event of (BM4) is inside F_k). By the Markov property at time n and at the stage ends, and Lemma 1 with
the remark of section 0, the probability in the statement is at least
q^J·(g − 2n·π(A_0^c)) − Σ_(k≤J) 2n·π(A_k^c),   g := P_π(cyl(t) ∩ S(n) ∩ {z_n ∈ box_0}) ≥ n^(−μ−ε/4).
q^J ≥ (δ^J)^((1+μ)(1+κ_δ)) ≥ n^(−(1+μ)−ε/4) since δ^J = r_J²/n ≥ 1/n. Choose M with ε_r²·M/4 ≥ 2μ + 10; then by Lemma 1
2n·π(A_k^c) ≤ n^(−2μ−5) for k ≤ J, and the claim follows for n large. ∎

**Corollary A1 (reduction to the end game).** With the parameters of Theorem A put
ρ(h) := inf e_N(s, z) over s ∈ A(h; ε_r, η₀, 1), z ∈ [ah, bh]² (floors as in section 0), N ∈ {⌊h²⌋, ⌊h²⌋ + 1}, N ≡ level(s).
Then R(2n) ≥ n^(−1−2μ−ε)·ρ(r_J)/π(cyl(t)), with M·log n ≤ r_J < M·log n/√δ. Consequently:

**(E) (end game; open).** *For every fixed choice of (δ, ε_r, η₀, a, b): ρ(h) ≥ exp(−o(h)) as h → ∞.*

(E) implies R(2n) ≥ n^(−α−o(1)), hence with Theorem U of 632, M_m = Λ^m·m^(−α+o(1)). The trivial bound is
ρ(h) ≥ exp(−C·h) (brute force, (F4)); that is Theorem Lo. A polynomial bound ρ(h) ≥ h^(−C) is what one expects.
(E) is weaker than statement (M) of 632 in three ways: only one scale is involved; no exponent has to be sharp
(any exp(−o(h)) will do); no classes below the scale log n have to be constructed. It is stronger in one way: it
asks for all states of A(h), whose top εh letters are arbitrary.

*Proof of the corollary.* Markov property at the end of stage J; the number of points left is r_J² up to the shift
by one, with the right parity (section 0). ∎

## 3. Two unconditional statements: the bottom of the words, and typical states

### 3.1 Words that are regular seen from the bottom

Letters of a word w of a true-chain state are indexed from the bottom: w_1 is the bottom letter, w_(abs(w)) the top.
For ε > 0, C_w ≥ 1, ρ₀ ≥ 1 the word w is *(ε, C_w, ρ₀)-regular* if for all 0 ≤ i < j ≤ abs(w) with j ≥ ρ₀ and
j − i ≥ C_w·log(2 + j):   abs(w_(i+1) + … + w_j) ≤ ε·(j − i).
A state is regular if both words are. 𝓡 = 𝓡(ε, C_w, ρ₀) is the class of regular states. The resolution is logarithmic
in the height above the bottom. The class is larger for larger C_w or ρ₀.

**(R1)** If w is regular and 0 ≤ b ≤ ℓ := abs(w), the sum of the top b letters is at most ε·b + C_w·log(2 + ℓ) + ρ₀ in
absolute value. (If ℓ < ρ₀ then b < ρ₀. If ℓ ≥ ρ₀ and b ≥ C_w·log(2 + ℓ), the window (ℓ − b, ℓ] is constrained. Otherwise
b < C_w·log(2 + ℓ).) So a regular word of length at most C·r is ε-regular at scale r in the sense of 632 as soon as
C_w·log(2 + Cr) + ρ₀ ≤ εr; in particular Lemma 2 (a) of 632 holds with ε replaced by 2ε on the event that the words stay
regular, for r ≥ r₀(ε, C_w, ρ₀).

**Proposition G (bottom regularity can be built into the growth; proved).** Fix ε_w > 0 and ε > 0. There are C′ and
0 < a < b such that for every C_w ≥ C′ there is ρ′ such that for every ρ₀ ≥ ρ′ and n ≥ n₀:

    P_π( cyl(t) ∩ S(n) ∩ {z at time n ∈ [a√n, b√n]²} ∩ {the two words above t at time n are (ε_w, C_w, ρ₀)-regular} ) ≥ n^(−μ−ε).

So statement (BR) of 632 is not needed as a hypothesis: the lower bound may choose the growth event, and the
growth event can include regularity of the lowest letters at every scale.

*Proof.* Follow 632 section 2.2 (scales r_j = √n·δ^((J−j)/2) increasing in j, times T_j = r_j², box_j = [a·r_j, b·r_j]²,
classes A_j, stage bound q_L ≥ cδ^μ/2) with two additions.
(i) *Freezing.* In the one-stage event require in addition that both heights stay ≥ (a/2)·r_j during stage j + 1
(from time T_j to T_(j+1)). In Brownian terms the cone Q is replaced by Q + (a/2)√δ(1,1) for a start at √δζ, ζ ∈ [a, b]²;
after translation this is (BM2) for starts in √δ·[a/2, b]², so the stage bound keeps the form c·δ^μ − ε′ with a smaller
c. Consequence: the letters at heights ≤ (a/2)·r_j above t do not change after time T_j.
(ii) *Regularity.* Let m_j := C_w·log(2 + (a/2)·r_(j−1)) and L_j := ⌈b·r_j⌉ + 1, and let Bad_j be the event that at some
move k ≤ T_j some window of at least m_j consecutive letters among the top L_j letters of σ or τ at move k has a
letter sum larger than ε_w times its length. At a fixed move the words have law π (a top letter and fair coins),
so by Hoeffding and a union bound P_π(Bad_j) ≤ 4·T_j·L_j²·exp(−(ε_w·m_j − 1)²/(2m_j)) ≤ C_δ·r_j^(4 − ε_w²·C_w/4).
Let Reg_j be the event that the words above t at time T_j are regular. Claim: Reg_(j−1), the freezing during stage j,
the box at time T_j and the complement of Bad_j imply Reg_j. Indeed a window (i, j′] of the word at time T_j with
j′ ≤ (a/2)·r_(j−1) lies in the frozen part and is constrained by Reg_(j−1); a window with j′ > (a/2)·r_(j−1) has length
≥ C_w·log(2 + j′) ≥ m_j and lies among the top L_j letters. Start: choose ρ₀ > b·r_(j₀), so that Reg_(j₀) is empty of
content (j₀ is the first stage of the induction of 632, r_(j₀) a constant depending on δ).
With g′_j := P_π(cyl(t) ∩ S(T_j) ∩ box_j ∩ Reg_j ∩ freezing so far), the recursion of 632 becomes
g′_(j+1) ≥ q_L·(g′_j − 2T_j·π(A_j^c)) − P_π(Bad_(j+1)), and g′_(j₀) = g_(j₀) ≥ exp(−c₁·r_(j₀)/4) as there. The induction
g′_j ≥ g′_(j₀)·(q_L/4)^(j−j₀) =: c·r_j^(−K) (K depends on δ) goes through if P_π(Bad_(j+1)) ≤ (q_L/4)·c·r_j^(−K), which holds for
all j ≥ j₀ when ε_w²·C_w/4 ≥ K + 5 and r_(j₀) is large. The end is as in 632: g′_J ≥ n^(−μ−o_δ(1)−o(1)). ∎

Remark. The same proof gives more: the frozen part of the words at each scale can be required to lie in any class
whose complement has stationary probability at most r^(−K−5) at scale r. The cost so far in the growth is only a
power of the current scale; this is why the growth is easy and the descent is not.

### 3.2 The end game from typical states

**Proposition T (proved).** For every ε > 0 there are c_ε > 0, C_ε, h_ε and classes B_h of bulk states with
π(B_h^c) ≤ C_ε·h^(C_ε)·exp(−c_ε·h), such that e_N(s, z) ≥ exp(−ε·h) for all h ≥ h_ε, s ∈ B_h, z ∈ [ah, bh]² (a floor of s)
and N ∈ {⌊h²⌋, ⌊h²⌋ + 1} with N ≡ level(s).

*Proof.* Use the constant-probability stages of Theorem Lo of 632 (δ = 1/16, stage probability q, classes
A(r) = A(r; ε_r, η₀, 1), last-scale bound: from a state of A(r) with floor in [ar, br]² the chain returns to t in
exactly the remaining ≈ r² points with probability ≥ q′·exp(−C₅·r)). Let ε₁ = ε/(2C₅), J′ the least integer with
δ^(J′/2) ≤ ε₁, r_i = δ^(i/2)·h. Define classes backward: B^(J′) := A(r_(J′)), and for i < J′

    B^(i) := { s ∈ A(r_i) : P_s( X_k ∉ B^(i+1) for some move k ≤ r_i² ) ≤ q/2 }.

For s ∈ B^(i) with floor in the box of scale r_i: P_s(stage i + 1 succeeds and ends in B^(i+1)) ≥ q − q/2. So for
s ∈ B_h := B^(0): e_N(s, z) ≥ (q/2)^(J′)·q′·exp(−C₅·ε₁·h) ≥ exp(−ε·h) for h large. By Markov's inequality and stationarity,
π(B^(i) c) ≤ π(A(r_i)^c) + (2/q)·r_i²·π(B^(i+1) c), so π(B_h^c) ≤ (J′ + 1)·(2h²/q)^(J′)·max_i π(A(r_i)^c), and by Lemma 1 of 632
π(A(r_i)^c) ≤ C·h³·exp(−ε_r²·δ^(1/2)·ε₁·h/4). ∎

So (E) holds for all states outside a set of exponentially small stationary measure, with rate c_ε of order ε.
**This does not combine with Theorem A.** To use B_h at the end of the descent in absolute probability one needs
c_ε·h ≥ (α + 1)·log n, that is h ≥ C·log n/ε, and then exp(−εh) = n^(−C′): the product ε·h cannot be made small. This is
the inequality of section 1 again, and it locates the gap exactly: *typical states are fine; what is missing is
that the state reached by the conditioned chain at height C·log n is typical at scales below its height, a
statement about a conditional law.*
By Proposition G (and its remark), the part of that state which was written during the growth can be required
typical in any such sense. The part that is not controlled is the layer of letters pushed during the last stages
of the descent and still present: of the order of log n letters.

## 4. Pushed letters: a comparison inequality, the statement (DR*), and the reduction

### 4.1 The chain from an arbitrary state against the stationary chain

For a sequence γ of k moves (each move is one entry of the move table) let a(γ), a′(γ) be the numbers of letters of
the start state popped along γ from σ and from τ; they are functions of γ alone. β = m₀/(1+2m₀).

**Proposition D (proved).** For every bulk state s, every k, d and every set Γ of sequences of k moves,

    P_s( Γ ∩ {a + a′ ≤ d} ) ≤ (4/β)·2^d·P_π( Γ ∩ {a + a′ ≤ d} ).

*Proof.* The probability that the bulk chain started at s performs γ depends only on the level of s and on its top
a(γ) + 1 and a′(γ) + 1 letters (625, proof of Corollary 1). Let cyl(γ, s) be the cylinder of states that agree with s
there. Then P_π(γ and X₀ ∈ cyl(γ, s)) = π(cyl(γ, s))·P_s(γ), and π(cyl(γ, s)) ≥ β·2^(−(a+a′+2)) by 625 Theorem 1. Sum over
γ ∈ Γ with a + a′ ≤ d. ∎

**Corollary D1 (dynamic regularity while the chain digs little; proved).** For every bulk state s and all k, ℓ, d, ε
with εℓ ≥ 4: the probability under P_s that at some move j ≤ k one has a_j + a′_j ≤ d and some window of at least ℓ
consecutive letters, among the letters pushed during the first j moves and present at move j, has a letter sum
larger than ε times its length, is at most (16/β)·k³·2^d·exp(−ε²ℓ/4).

*Proof.* The event is a set of move sequences (pushed letters are functions of the moves). Under π the word at a
fixed move is a top letter followed by fair coins; Hoeffding, a union bound over at most 4j² windows and over
j ≤ k, then Proposition D. ∎

So (DR) of 632 holds from every state, with no hypothesis on its letters, as long as the number d of revealed
letters of the start state is at most ε²ℓ/8 − C·log k. It does not give (DR): in a stage at scale r the chain reveals
of the order of r letters and the windows that matter have length C·log r. The factor 2^d is the price of
conditioning on the exact revealed letters; a proof of (DR*) below has to use the revealed letters only through
their statistics.

### 4.2 The statement (DR*)

**(DR*) (conjecture).** For every ε > 0, η > 0, C′, ρ′ there are C_w ≥ C′ and ρ₀ ≥ ρ′ such that for every r ≥ 1 and every
true-chain state s ∈ 𝓡(ε, C_w, ρ₀) with both heights at most r:

    P_s( for some move k ≤ r² the chain is alive after k moves and X_k ∉ 𝓡(ε, C_w, ρ₀) ) ≤ η.

In words: from a state that is regular seen from the bottom at logarithmic resolution, the chain stays in that
class during r² moves, up to a small conditional probability, uniformly in r and in the state. Under π the
corresponding absolute statement is the computation in the proof of Proposition G. Evidence: section 5.3 (data).
(DR*) is a quenched statement: the buried letters of s are fixed, and only their window sums are known.

### 4.3 The invariance principle from regular states, for the killed chain

**Lemma 2′.** Fix u ∈ [1/2, 1], a closed box D ⊂ Q, θ > 0 and ε′ > 0; D^(−θ) is the set of points of D at distance ≥ θ
from its boundary. There are ε > 0, η > 0 such that for all C_w, ρ₀ there is r₀ with the following property. Let r ≥ r₀,
let s ∈ 𝓡(ε, C_w, ρ₀) be a true-chain state with heights z⁰ ∈ [0, br]², and suppose that the event G′ := {X_k ∈ 𝓡 at every
move k ≤ r² after which the chain is alive} has P_s(G′ fails) ≤ η. Then

    P_s( alive during u·r² points, heights at the end in r·D ) ≥ P_(z⁰/r)( W ∈ Q + θ(1,1) on [0, u], W_u ∈ D^(−θ) ) − ε′ − η.

*Proof.* As Lemma 2 (c) of 632, by contradiction, with one new device for the killing. Suppose the statement fails:
there are sequences ε_n, η_n → 0, constants C_(w,n), ρ_(0,n), scales r_n with (C_(w,n)·log r_n + ρ_(0,n))/r_n → 0 and states
s_n violating the inequality (write u_n = u, D_n = D, θ_n = θ below). Realize the true chain from s_n as the bulk chain on the words
of s_n placed on top of the letters of t, and let ζ be the first move at which a letter of t is popped (killing).
Before ζ, (F3) reads L_k = J·z_k + (1/2)δ^σ·D^σ_k + (1/2)δ^τ·D^τ_k + O(1) with D^σ_k = S(σ_k; b_k) − S(σ_0; a_k), where the popped
letters are letters of the word of s_n. On G′, by (R1), abs(D^σ_k) ≤ ε_n·(a_k + b_k) + 2C_(w,n)·log(2 + C·r_n + C·Z*_k) + 2ρ_(0,n),
the same for τ; so Lemma 2 (a) of 632 holds up to ζ with an error that is o(r_n)·(1 + L*/r_n).
On an enlarged probability space take a Brownian motion B′ with covariance H independent of the chain and put
M̃ⁿ(v) := r_n⁻¹·L_(k(r_n²v)∧ζ) + B′(v) − B′(v ∧ ζ_n), ζ_n := T_ζ/r_n². M̃ⁿ is a martingale with jumps ≤ 3/r_n; its predictable
quadratic variation is r_n⁻²·A_(k(r_n²v)∧ζ) + H·(v − v ∧ ζ_n). By Lemma 0 of 632 the first term is
H·(v ∧ ζ_n) + r_n⁻²·(c·N + γ·Δ + R) + O(r_n⁻²) at the stopped time; N is bounded by C·(r_n + Z*) on G′, Z*/r_n and
max abs(R)/r_n are tight (Doob for the stopped martingales, at most r_n² moves). So the quadratic variation tends to
H·v in probability and, by (O2), M̃ⁿ converges in law to a Brownian motion with covariance H. Put
Ỹⁿ := z⁰_n/r_n + J⁻¹·M̃ⁿ. Along a subsequence z⁰_n/r_n → ξ, and Ỹⁿ converges in law to ξ + W.
On G′, sup over v ≤ min(u_n, ζ_n) of abs(r_n⁻¹·(z⁰_n + z_(k(r_n²v))) − Ỹⁿ(v)) → 0 in probability (including the killing move,
which changes a height by one). If Ỹⁿ stays in Q + (θ_n/2)(1,1) on [0, u_n] and this supremum is below θ_n/4, then
ζ_n > u_n: at the killing move a true height is negative, while Ỹⁿ is at least θ_n/2 there. So
P_(s_n)(alive during u_n·r_n² points, end in r_n·D_n) ≥ P(Ỹⁿ ∈ interior of Q + (θ_n/2)(1,1) on [0, u_n], Ỹⁿ(u_n) ∈ interior of
D_n^(−θ_n/2)) − η_n − o(1). By the portmanteau theorem for open sets, and because the open event with θ/2 for the start ξ contains the
closed event with θ for every start close to ξ, the right side is at least P_(z⁰_n/r_n)(W ∈ Q + θ_n(1,1) on [0, u_n], W_(u_n) ∈ D_n^(−θ_n)) − η_n − o(1). Contradiction. ∎

### 4.4 The reduction

**Theorem B (reduction; proved).** (DR*) implies R(2n) ≥ n^(−α−o(1)), hence, with Theorem U of 632,
M_m = Λ^m·m^(−α+o(1)) for the band [−1, 2].

*Proof.* Fix ε > 0. Take δ, θ = a√δ/8 and c₄ as in Theorem A, with κ_δ ≤ ε/(4(1+μ)); ε′ := c₄δ^(1+μ)/8 (so ε′ ≤ θ); (ε_reg, η)
from Lemma 2′ for (u, D, θ, ε′) = (1 − δ, √δ·[a, b]², θ, ε′), with η ≤ c₄δ^(1+μ)/8; (C′, ρ′) from Proposition G for ε_w = ε_reg and exponent ε/4; (C_w, ρ₀) from (DR*)
for (ε_reg, η, C′, ρ′); r₀ from Lemma 2′. 𝓡 := 𝓡(ε_reg, C_w, ρ₀). Stages and boxes as in Theorem A, now for all k ≤ J with J
the largest k such that r_k ≥ r₀.
*One stage.* Let s ∈ 𝓡 have heights in box_(k−1), r_(k−1) ≥ r₀. By (DR*) with r = b·r_(k−1), P_s(G′ fails) ≤ η. By Lemma 2′
(scale r_(k−1), u = 1 − δ, D = √δ·[a, b]², the θ and ε′ just fixed; r₀ is taken ≥ 5b/√δ) and (BM4), P_s(F_k) ≥ c₄δ^(1+μ) − ε′ − η, hence
P_s(F_k and the state at the end of stage k is in 𝓡) ≥ c₄δ^(1+μ) − ε′ − 2η ≥ q := c₄δ^(1+μ)/2.
*All stages.* By the Markov property, for every s ∈ 𝓡 with heights in box_0: P_s(F_1 ∩ … ∩ F_J) ≥ q^J ≥ n^(−(1+μ)−ε/4).
*The end.* After stage J the heights are at most b·r₀/√δ and r_J² ≤ r₀²/δ points remain, with the right parity. By
(F4) the chain reaches t in at most 4b·r₀/√δ + 4 points and fills the rest with pairs of gos "nothing changes" at t;
this has probability at least a constant c(δ) > 0. So e_N(s, z) ≥ c(δ)·n^(−(1+μ)−ε/4) for the true-chain states s ∈ 𝓡 with
heights z ∈ box_0 and the right N ∈ {n, n + 1}.
*The growth.* By Proposition G the event that the chain started at t survives n points, has heights in box_0 and
is in 𝓡 at time n has P_π(· ∩ cyl(t)) ≥ n^(−μ−ε/4). The Markov property at time n gives
R(2n) = P_π(cyl(t) ∩ E_(2n))/π(cyl(t)) ≥ c(δ)·n^(−1−2μ−ε/2)/π(cyl(t)) for n large. ε is arbitrary. ∎

Compared with 632 Proposition R: statement (BR) is gone (Proposition G), the classes are explicit, and the
exceptional probability in the dynamic statement only has to be a small constant η, not exp(−c·ℓ): a failure of
regularity is charged to the stage in which it happens.
What Theorem B uses: 632 Lemma 0, (F1) to (F4), (O1), (O2), (BM2), Theorem S (through Proposition G), Theorem U;
625 Theorem 1. It does not use Lemma 1 of 632 except inside Proposition G.

## 5. The driving symbols, exact flux identities, and what the other inputs give

### 5.1 The conditional law of the pushed letters (idea 1 of the brief)

The bulk chain is a deterministic function of the start state and of an i.i.d. sequence of symbols
(class, three coins): class A with probability 1/(1+m₀) or B; a coin "lose" (probability 2m₀/(1+2m₀)), a coin
"receive −" (m₀/(1+m₀)), a coin "receive +" (1/(1+2m₀)). At active top +: A is a stay (the active stack loses its top if
"lose", the passive stack receives − if "receive −"), B is a go (the passive stack receives + if "receive +"). At
active top −: every symbol is a go; A pops the top; the passive stack receives + if "receive +".
So a pushed letter is − exactly when it is pushed by a stay, and whether a symbol A is a stay is decided by the
active top, which is a letter pushed earlier or a letter of the start state. Given the start state the sequence
of pushed letters is therefore not a hidden-Markov sequence with a finite hidden state: the hidden state is the
whole pair of stacks. Two exact statements survive.

(i) *Uniform ellipticity.* Conditionally on the whole past, at every move at level 1 a letter + is pushed onto τ
with probability ≥ (m₀/(1+m₀))/(1+2m₀), and within three moves a letter − is pushed onto τ before any + with
probability bounded below; the same for σ. So no pattern of j pushed letters has conditional probability above
(1 − c)^j. This is far from balance and gives no large-deviation bound for letter sums.

(ii) **Lemma F (flux identities; proved).** For the bulk chain from any state let G^(1+)_k, G^(1−)_k, G^(0+)_k, G^(0−)_k be the
numbers of gos among the first k moves made at level 1 or 0 with active top + or −, G¹ = G^(1+) + G^(1−), G⁰ likewise,
G^± = G^(1±) + G^(0±). Then the following are martingales with bounded increments, from every start state:

    N_(σ+) − (G⁰ − 2G^(1+))/(1+2m₀),      N_(σ−) − (G^(0+) − G^(1−))/(1+m₀),
    N_(τ+) − (G¹ − 2G^(0+))/(1+2m₀),      N_(τ−) − (G^(1+) − G^(0−))/(1+m₀),
    x + y − c₀·(G^+ − G^−),      c₀ = m₀/((1+m₀)(1+2m₀)) = 0.16877….

*Proof.* Conditional drifts in the four states (level, active top), from the table of 630 section 4.2. Level 1,
top +: a stay has probability 1/(1+m₀), a go m₀/(1+m₀); E[ΔN_(σ+)] = −(1/(1+m₀))·2m₀/(1+2m₀) and E[ΔG^(1+)] = m₀/(1+m₀), so the
first expression has drift −2m₀/((1+m₀)(1+2m₀)) + 2m₀/((1+m₀)(1+2m₀)) = 0; E[ΔN_(τ−)] = (1/(1+m₀))·m₀/(1+m₀) = E[ΔG^(1+)]/(1+m₀);
E[ΔN_(τ+)] = (m₀/(1+m₀))/(1+2m₀) = E[ΔG¹]/(1+2m₀). Level 1, top −: one go; E[ΔN_(σ−)] = −1/(1+m₀) = −E[ΔG^(1−)]/(1+m₀);
E[ΔN_(τ+)] = 1/(1+2m₀) = E[ΔG¹]/(1+2m₀); the other two expressions do not move. Level 0 is the mirror image. The last
line is the sum of the four. `flux_check.py` confirms all twenty drifts to 40 digits. ∎

Consequences. (a) The total height is, up to a martingale, c₀ times (gos made at top +) − (gos made at top −): tops +
make the stacks grow, tops − make them shrink. (b) The two martingales of 625 are the combinations in which the
go counts cancel (G¹ − G⁰ is the level change). (c) The letter imbalances D^σ, D^τ of 625 Lemma 3 are, up to martingales,
linear in the two numbers G^(1+) − G^(1−) and G^(0+) − G^(0−). So **(DR*) is equivalent to a statement about the top letter
only**: from a regular state, in every time window the gos at each level are made at top + and at top − equally
often, up to ε·(window sum scale) and the martingale terms. Under π this is the identity P(go at top +) = P(go at
top −) = m₀/(1+2m₀) per move. I did not find a proof from a fixed state.

### 5.2 The chain without killing of 633 (Corollary 3.4), and time reversal

*633.* P*(s, s′) = P_c(s, s′)·h(s′)/h(s) with h = V*/V_c, h(t) = 1; h depends on the two bottom runs only (factor
2^(r+1) − 1 for a bottom run of σ of length r, and the analogous factor for τ). Hence P*_s(at t after N points) =
e_N(s)/h(s): the return probabilities to t are the same numbers, and a bound for one chain is a bound for the
other. Two remarks. (1) h is not the harmonic function of the cone: it is bounded on states whose bottom runs
are bounded and grows exponentially in the bottom runs, while the function that describes conditioning on long
survival grows like a power 2μ of the heights. So P* is not the chain conditioned to survive, and it gives no
information on the law of the lowest letters under that conditioning. (For a long bottom run of σ a push of + onto
it is favored by a factor tending to 2 and a pop is disfavored by 2; presumably P* escapes by growing that run
linearly. Not analyzed.) (2) In the lower bound the lowest letters are handled by Proposition G, and the open
point (DR*) concerns the top layer, where P* and P_c agree. So Corollary 3.4 of 633 does not shorten the present
argument. It is not used.

*Time reversal (idea 4).* Not derived. The Markov property at time n gives the pairing
R(2n)·π(cyl(t)) = Σ_s P_π(cyl(t), S(n), X_n = s)·e_n(s, z(s)), which is not a sum of squares. I checked on the
move table that the time reversal of the stationary bulk chain is not the same macro chain after exchanging the
stacks and flipping the letters: in the reversed chain a move is decided by the top of the stack that was
passive (a top − can only be un-pushed by reversing a stay), and for moves in which the active stack kept its
top also by that top. This agrees with 633 (no bijection of macro states). A reversal at the level of raw
machine states with the river edge exists by the symmetry of the band, and would turn e_n(s, z) into a growth
probability; but a lower bound for the pairing still needs the two factors to be large on a common set of words,
which is the same quenched question. Status of this paragraph: remark.

### 5.3 Data for (DR*) (`dr_probe.py`, `dr_probe.log`; not a proof, not a fit)

Bulk chain from deterministic bottomless stacks with a buried pattern, 50,000 moves, 400 runs per pattern. S is the
letter sum and b the number of the letters pushed during the run and present at its end; the table gives
(sum of S over the runs)/(sum of b over the runs) with its standard error, for σ and for τ.

| buried pattern | σ | τ | mean b (σ, τ) |
|---|---|---|---|
| fair coins (control) | −0.001 ± 0.004 | 0.003 ± 0.004 | 124, 131 |
| (+−) periodic | 0.003 ± 0.005 | 0.006 ± 0.005 | 112, 123 |
| (−+) periodic | −0.002 ± 0.004 | −0.001 ± 0.005 | 123, 122 |
| (++−−) | 0.005 ± 0.004 | −0.004 ± 0.005 | 123, 126 |
| (+++−−−) | 0.004 ± 0.005 | −0.005 ± 0.005 | 120, 125 |
| (++++++−−−−−−) | 0.004 ± 0.005 | −0.002 ± 0.005 | 122, 117 |
| all + | −0.029 ± 0.005 | −0.040 ± 0.005 | 157, 169 |
| all − | −0.001 ± 0.005 | −0.005 ± 0.005 | 99, 93 |

For every buried pattern that is balanced in windows the pushed layer is balanced within the error (0.5 percent),
whatever the order of the letters. This is consistent with (DR*). A buried stack of equal letters + (which is not
in the class 𝓡) gives a pushed layer with an excess of − of 3 to 4 percent at this size; the run does not decide
whether that is a boundary effect of a few letters (S of order −5) or a proportional one. The probe tests means
only, not the small probabilities in (DR*).
