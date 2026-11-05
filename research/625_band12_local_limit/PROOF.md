# 625: the band [−1, 2]. Stationary law, exact growth rate, invariance principle, and the lemma that is still missing (2026-10-10)

Status words as in `WAVE_620_COMMON.md`: **proved** (complete argument here), **sketch**, **conjecture**, **data**.

**Inputs that are not proved here.**
- AE1 (label machine, `606_rung2_feasibility/NOTES.md` section 1; rests on Lemmas 1 to 3 of `204_winding_bridge/SKETCH.md`;
  unrefereed).
- Propositions 2 and 3 of `622_rung2_stack_method/PROOF_BAND12.md` (ledger AH2, AH3; proved there, not reviewed).
  I did not re-derive Proposition 2 from the machine rules. My own program for the move table (`projtest.py`)
  reproduces M_m·z_c^(2m−2) exactly for m ≤ 10, and the four equations of Proposition 3 at the critical member are
  rechecked exactly (`qm.py`).

Every statement about M_m below is conditional on these. Statements about the chain P_c itself (sections 2, 3, 5)
are unconditional: they are statements about an explicitly defined Markov chain.

## 0. Setting and notation

M_m = M^[−1,2]_m. m₀ = 0.917543340819817684… is the positive root of 4x⁴ − 2x − 1. All exact computations are done in
the field Q(m₀) (`qm.py`).

Critical weights: p = p_+ = q_+ = 1/(2m₀), p_− = q_− = m₀, κ = 1, z = z_c = 1/(2m₀²(1+m₀)),
Λ = 1/z² = (1+2m₀)(1+m₀)² = 10.4245356612….

State: (level, σ, τ), σ and τ words over {+, −}, written with the top first: σ = s₁s₂…, τ = t₁t₂…. At level 1 the stack
σ is *active* and τ is *passive*; at level 0 the roles are exchanged. One *move* of the critical chain P_c, at a state
whose active stack has top letter a (Proposition 2 of 622 with the weights z^(points)·V(s′)/V(s)):

    a = + :   stay (2 points):  active keeps / loses its top, passive keeps / receives − :
                                z², z²m₀, z²/p, z²m₀/p
              go (1 point, level changes): active unchanged, passive keeps / receives + :      z, zp
    a = − :   go (1 point): active keeps / loses its top, passive keeps / receives + :         z, zp, z/m₀, zp/m₀

The probabilities sum to 1 for both a (exact, `qm.py`). Useful values: z(1+p) = m₀/(1+m₀), so

    a = + :  P(stay) = 1/(1+m₀),  P(go) = m₀/(1+m₀);        a = − :  P(pop) = 1/(1+m₀),  P(no pop) = m₀/(1+m₀).

Three chains are used.
- The **true chain**: the machine started at t = (level 0, σ empty, τ = (+)), killed when a pop is not available (σ
  empty at level 1 allows only the two moves "go, σ unchanged"; the bottom letter of τ is never removed).
  Exact identity (622, (4.1)): M_m = Λ^(m−1) · P_c(the true chain started at t is at t after 2m − 2 points).
- The **bulk chain**: the same move table on bottomless stacks, state space Ω = {0,1} × {+,−}^N × {+,−}^N. It is a
  Markov chain with a probability kernel. At a state of the true chain whose active stack has a removable top, the
  true kernel and the bulk kernel agree.
- The **stationary bulk chain**: the bulk chain started from the law π of Theorem 1.

Time is counted in moves (index k) or in points: T_k = number of points of the first k moves, k ≤ T_k ≤ 2k.
x_k, y_k are the heights of σ, τ relative to time 0 (pushes minus pops), well defined for the bulk chain.
N_c(k), c ∈ {σ+, σ−, τ+, τ−}, is the number of letters of kind c pushed minus the number popped during the first k
moves; x = N_(σ+) + N_(σ−), y = N_(τ+) + N_(τ−).

## 1. The hidden-letter lemma in its pathwise form is false [proved; exact computation]

Statement tested: "conditionally on the history of (level, heights, two top letters), a letter revealed by a pop is a
fresh fair coin, so the projection to (level, heights, tops) is a Markov chain".

**It fails, for a structural reason.** A letter records the kind of move that pushed it: at level 1 a stay pushes −
onto τ and a go pushes +. The history of (level, heights) therefore determines every letter in both stacks. A letter
revealed by a pop is a function of the observed history, not a fresh draw. Two histories with the same current
projected state can have different second letters, hence different laws for the next projected state. So the
projection to (level, heights, tops) is not a Markov chain.

**The weaker hope fails too** (that the fresh-letter chain has the same return probabilities). `projtest.py` computes,
in exact arithmetic in Q(m₀), the return probability to t of the projected chain with fresh letters (+ with
probability b) and compares it with M_m·z^(2m−2):

| m | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 |
|---|---|---|---|---|---|---|---|---|---|---|
| projected / true, b = 1/2 | 1 | 1 | 1 | 0.993450 | 0.991446 | 0.995777 | 1.005074 | 1.017594 | 1.031973 | 1.047248 |

Equality holds exactly for m ≤ 3 and fails exactly from m = 4 on. No other coin repairs it: for b = 1/4, 2/5, 9/20,
11/20, 3/5, 3/4 the ratios at m = 4..10 also differ from 1 and drift (`projtest.log`). The drift of the ratio suggests
that the fresh-letter walk, which is a genuine finite-mode quarter-plane walk, has a different exponent (data, seven
terms). It is not the band.

The correct statement is about the stationary law at a fixed time, not about the conditional law given the history.

## 2. The stationary law of the bulk chain [proved]

**Theorem 1.** Let π be the probability measure on Ω under which
- the level is uniform on {0, 1};
- all letters of both stacks are independent of each other and of the level;
- every letter is a fair coin, except the top letter of the active stack, which is − with probability m₀/(1+2m₀) =
  0.323638… and + with probability (1+m₀)/(1+2m₀).

Then π is invariant for the bulk chain.

*Proof.* Let X₀ have law π and let X₁ be the state after one move. The move is a function of the level, the active
top and an independent uniform variable. Each stack changes by at most one letter at its top.

*Step 1 (reduction to a finite check).* Write W₁ for the vector (level, s₁, s₂, t₁, t₂) at time 1. Fix the move and
the values of (level, s₁, t₁) at time 0. If σ keeps its top, the letters of σ at depth ≥ 3 at time 1 are the letters at
depth ≥ 3 at time 0. If σ loses its top, they are the old letters at depth ≥ 4. If σ receives a letter, they are the old
letters at depth ≥ 2. In each case they are old letters at depth ≥ 2 that do not enter W₁ (W₁ uses the old letters at
depth ≤ 2, ≤ 3, ≤ 1 respectively). Under π the old letters at depth ≥ 2 are independent fair coins, independent of the
move and of the old tops. The same holds for τ. Hence, at time 1, the letters at depth ≥ 3 of both stacks are
independent fair coins and independent of W₁. So the law of X₁ is (law of W₁) × (fair coins), and π is invariant if
and only if W₁ has the law that π gives to (level, s₁, s₂, t₁, t₂).

*Step 2 (the finite check).* The law of W₁ is computed from the law that π gives to (level, s₁, s₂, s₃, t₁, t₂, t₃) and
the move table: 128 starting windows, at most 6 moves each, 32 target cylinders. `stationary.py` does this in exact
arithmetic in Q(m₀) and finds equality on all 32 cylinders. ∎

The law of the active top can be found by hand: with u, v the stationary probabilities of (level 1, top +) and
(level 1, top −), the rate of entering level 1 is G = u·m₀/(1+m₀) + v, a new active top after a go is the pushed + with
probability p/(1+p) and an old passive (fair) letter otherwise, and after a pop it is a fair letter. The balance
equation for v reads 1/2 − u = w/4 with w = 1/(1+p) = 2m₀/(1+2m₀), which gives v = m₀/(2(1+2m₀)).

Remarks.
- This closes the gap named in 622 section 4 (ii): the buried letters are independent fair coins, and the "factor
  depending on the tops" is a bias of the active top only.
- Mean number of points per move under π: 1 + P(stay) = 1 + 1/(1+2m₀) = (2+2m₀)/(1+2m₀) = 1.352722….
- Net stationary flux of each of the four kinds of letters is zero (exact, `stationary.py`); in particular both heights
  have zero drift.

**Corollary 1 (fresh letters at new minima).** For the stationary bulk chain let a_k = −min_(j≤k) x_j be the number of
time-0 letters of σ popped during the first k moves, and a′_k the same for τ. Conditionally on the first k moves, the
time-0 letters of σ at depth > a_k + 1 and of τ at depth > a′_k + 1 are independent fair coins.

*Proof.* By induction on k, the event that the first k moves are a given sequence is determined by the level at time
0, the uniform variables, and the time-0 letters at depth ≤ a_k + 1 and ≤ a′_k + 1 (a move needs the current active
top, which is either a letter pushed earlier, hence a function of earlier moves, or a time-0 letter at depth at most
one more than the number of pops so far). Under π the deeper letters are independent of these. ∎

So the correct form of the hidden-letter lemma is: *a letter of the initial stack revealed for the first time is a
fresh fair coin; a letter pushed during the observation is not.* Consequences:
- At the times k with x_k = min_(j≤k) x_j and y_k = min_(j≤k) y_j (simultaneous running minima) the future is
  independent of the past given (level, the two top letters). These times form a Markov renewal set R with finitely
  many modes. This is the only regeneration structure I found.
- R does not give cycles in the sense of Theorem 2 of 227: by Theorem 3 below, the probability that there is no point
  of R in [εn, n] is bounded below for every ε > 0, so the gaps have no exponential moment (hypothesis (A1) fails), and
  in the Brownian limit simultaneous running minima should not occur at positive times, because the cone angle is
  below 90 degrees (the criterion for cone points of planar Brownian motion; cited from memory, not used in any proof).
- The excursion t → t never reveals an initial letter. Every letter it pops was pushed during the excursion.
  Corollary 1 therefore gives no Markov property inside the excursion.

**Exact identity with the stationary chain.** In the bulk chain let E_n be the event: the time-0 top letters of σ and τ
are never popped, and after n points the level is 0 and both heights are 0. Then

    M_m · Λ^(−(m−1)) = P_π(E_(2m−2) | level 0 at time 0, s₁ = −, t₁ = +).                                        (2.1)

*Proof.* On E_n the letters below the initial tops are never used. The initial top + of τ plays the role of the bottom
letter of τ in the true chain. The initial top − of σ plays the role of "σ empty": at level 1 it allows the two gos
without pop and the two gos with pop, and the latter are excluded by E_n. So the paths of E_n from this start are in
bijection with the paths t → t of the true chain, with the same probabilities. ∎

## 3. Two martingales, the letter imbalance, and a fluctuation bound [proved]

Parametrize the product family of 622 (Proposition 3) by u = (u₁, u₂) = (log p, log p_−): p_+ = q_+ = p,
q_− = p(1+p)/p_−, z² = p/((1+p)(1+p_−)(1+q_−)), κ² = p_−(1+q_−)/(q_−(1+p_−)). For every member the ratio
V_u(X_k)·z(u)^(T_k) / (V_c(X_k)·z_c^(T_k)) is a martingale of the bulk chain (it is the ratio of two Doob transforms;
622 section 4 (i)). Put e₁ = (1+2p)/(1+p) = 2(1+m₀)/(1+2m₀) at the critical member. The first derivatives of log z vanish
there, and the first derivatives of the martingale are

    L¹_k = N_(σ+) + N_(τ+) + e₁·N_(τ−) − (e₁/(2(1+m₀)))·Δ_k,
    L²_k = N_(σ−) − N_(τ−) + (1/(1+m₀))·Δ_k,            Δ_k = [level 0 at time k] − [level 0 at time 0].

**Proposition 2.** (a) L¹ and L² are martingales of the bulk chain, from every starting state, with increments
bounded by a constant B.
(b) (second-order identity) With H = −(Hessian of log z in u) at the critical member, for every state and i, j ∈ {1, 2},

    E[ΔLⁱ ΔLʲ | state] = H_ij · E[ΔT | state] − E[β_ij · ΔN + γ_ij · ΔΔ | state],

where β_ij = ∂_i∂_j (log p_+, log p_−, log q_+, log q_−) and γ_ij = ∂_i∂_j log κ are constants (β₁₁ = (0, 0, 0, p/(1+p)²),
β₁₂ = β₂₂ = 0).
(c) Under π: E[ΔLⁱ ΔLʲ] = t̄·H_ij per move, with t̄ = (2+2m₀)/(1+2m₀).

*Proof.* The conditional law of a move depends only on (level, active top): four cases. (a) and (b) are the first
and second derivatives in u of the identity "the probabilities of the member u sum to 1"; `stationary.py` verifies
both in all four cases in exact arithmetic (it computes the left sides from the move table and the right sides from
the closed forms of the derivatives). (c) follows from (b), Theorem 1 and the zero fluxes. ∎

Values: H₂₂ = m₀/(1+m₀)², H₁₂ = −e₁m₀/(2(1+m₀)²), H₁₁ = (1/2)[p/(1+p)² + e₁²m₀/(1+m₀)² + (m₀/(1+m₀))·p/(1+p)²].

**Lemma 3 (letter imbalance; exact identity).** For the bulk chain let b_k = x_k − min_(j≤k) x_j be the number of
letters of σ at time k that were pushed after time 0, and a_k as in Corollary 1. Write S(w; j) for the sum of the top j
letters of the word w (+ counts 1, − counts −1). Then

    D^σ_k := N_(σ+)(k) − N_(σ−)(k) = S(σ at time k; b_k) − S(σ at time 0; a_k),

and the same for τ.

*Proof.* A letter pushed and popped during the first k moves cancels. What remains is the letters pushed after time 0
and still present (the top b_k letters at time k) minus the time-0 letters that were popped (the top a_k letters at
time 0). ∎

Since N_(σ±) = (x ± D^σ)/2 and N_(τ±) = (y ± D^τ)/2,

    (L¹_k, L²_k)ᵀ = J·(x_k, y_k)ᵀ + (1/2)·δ^σ·D^σ_k + (1/2)·δ^τ·D^τ_k + (bounded term),
    J = [[1/2, (1+e₁)/2], [1/2, −1/2]],   δ^σ = (1, −1)ᵀ,   δ^τ = (1 − e₁, 1)ᵀ,   det J = −(2+e₁)/4 ≠ 0.               (3.1)

**Lemma 4 (fluctuation bound).** There are constants K and n₀ such that for the stationary bulk chain and all n ≥ n₀,

    P_π( max over k ≤ n of max(|x_k|, |y_k|) ≤ K·√n ) ≥ 1/2.

Moreover, for every ε > 0 there is K_ε with probability ≥ 1 − ε, and on an event of probability 1 − o(1) one has in
addition |D^σ_k|, |D^τ_k| ≤ 2n^(0.4) + 2 for all k ≤ n.

*Proof.* Let J_n = ⌈2K√n⌉ + 2 and let G be the event that |S(σ at time k; j)| ≤ n^(0.4) + 1 and |S(τ at time k; j)| ≤
n^(0.4) + 1 for all k ≤ n and j ≤ J_n. By Theorem 1 the stack at a fixed time k has law π: the top letter is bounded
by 1 and the others are fair coins. By Hoeffding's inequality and a union bound,
P(G^c) ≤ 8(n+1)·J_n·exp(−n^(0.8)/(2J_n)) → 0.
By Proposition 2 (c) and Doob's inequality, P(max_(k≤n) |Lⁱ_k| ≥ K′√n) ≤ t̄·H_ii/K′².
Let η be the first k at which max(|x_k|, |y_k|) > K√n. If η ≤ n, then for k ≤ η we have a_k, b_k ≤ 2K√n + 2 ≤ J_n, so on G
Lemma 3 gives |D^σ_k|, |D^τ_k| ≤ 2n^(0.4) + 2, and (3.1) gives
|(x_η, y_η)| ≤ ‖J⁻¹‖·(√2·max_i max_k |Lⁱ_k| + C·n^(0.4) + C). On the event max |Lⁱ| < K′√n this is below K√n for
K = 2‖J⁻¹‖K′ and n large, a contradiction. So P(η ≤ n) ≤ 2t̄·max(H₁₁, H₂₂)/K′² + o(1), which is below 1/2 (or below ε) for K′
large. ∎

## 4. The growth rate [proved, given the inputs]

**Theorem 2.** There is a constant C such that for every m ≥ 1

    Λ^(m−1) · exp(−C·√m) ≤ M^[−1,2]_m ≤ Λ^(m−1),        Λ = (1+2m₀)(1+m₀)², 4m₀⁴ = 2m₀ + 1.

In particular lim M_m^(1/m) = Λ = 10.4245356612…, the root of 16x⁴ − 112x³ − 520x² − 536x + 25.

*Proof.* The upper bound is Theorem A of 622. For the lower bound it is enough, by M_m = Λ^(m−1)·P_c(t → t in 2m − 2
points) and M_m ≥ 1, to show P_c(t → t in n points) ≥ exp(−C√n) for all large even n. Let δ > 0 be the smallest
probability in the move table.

*(i) Any state can be built and removed at linear cost.* From t, with τ = (+): a stay with "σ receives −" appends − to
σ; a go with "σ receives +" followed by the go "nothing changes" (allowed at level 1 for every top) appends + to σ. So
every σ is built with τ = (+). Then, with σ fixed: a go from level 1 with "τ receives +" followed by the go "nothing
changes" from level 0 appends + to τ. To append − to τ, if the top of σ is +, use a stay with "σ keeps, τ receives −";
if not, first append a + to σ, use the stay, and remove the + by a stay with "σ loses its top". So every state
(level, σ, τ) with bottom letter + in τ is reached from t in at most C₁(|σ| + |τ| + 1) points, every move having
probability ≥ δ. Conversely, a top + of the active stack is removed by a stay with pop, a top − by a go with pop
followed by the go "nothing changes"; so t is reached from every state in at most C₁(|σ| + |τ| + 1) points.

*(ii) A good window.* Let K, n₀ be as in Lemma 4, T ≥ n₀, L = ⌈K√T⌉ + 2, and let F be the event that
max(|x_k|, |y_k|) ≤ L − 2 for all k ≤ T (T moves cover at least T points). Whether F occurs is determined by the level at
time 0, the top L letters of both stacks at time 0 and the uniform variables. P_π(F) ≥ 1/2, so there is a level and a
pair of words (w, w′) of length L with P(F | level, top windows w, w′) ≥ 1/2. Let s* = (that level, σ = w, τ = w′
followed by the bottom letter +). On F the bulk chain started with these windows never comes within two letters of
the bottom of the windows, so every state it visits has a removable active top, where the true kernel equals the bulk
kernel. Hence the true chain started at s* follows F with probability ≥ 1/2, and stops, at the first move at which
the number of points is ≥ T, in a state with both heights ≤ 2L + 1.

*(iii) Assembly.* Given a large even n, put T = n − 8C₁(K√n + 4) − 2. Go from t to s* (at most C₁(2L + 2) points), run
the chain as in (ii) (T or T + 1 points), go back to t (at most C₁(4L + 4) points), and fill the remaining number of
points, which is even (a path t → t has an even number of gos) and at most 8C₁(K√n + 4) + 2, with pairs of gos "nothing
changes" at t (probability z² per pair). The probability of the set of paths so described is at least
δ^(C₁(6L+6)) · (1/2) · z^(8C₁(K√n+4)+2) ≥ exp(−C√n). ∎

Remarks.
- 622 had the upper bound; the lower bound and hence the exact growth rate are new. Claim B of 622 is now a theorem.
- M_m·Λ^(−m) → 0: by (2.1) and Theorem 3 below in its mixing form, P_π(E_n | cylinder) ≤ P_π(x_k ≥ 0 for all k ≤ n/2 |
  cylinder) → 0. [proved given the mixing form of the martingale central limit theorem, see section 5.]

## 5. The invariance principle and the covariance [proved]

**Theorem 3.** For the stationary bulk chain let k(s) be the first k with T_k ≥ s (the number of moves needed for s
points). As n → ∞ the processes

    t ↦ n^(−1/2) · (x_(k(nt)), y_(k(nt))),    0 ≤ t ≤ 1,

converge in distribution (Skorokhod space) to a Brownian motion with covariance matrix Σ per point,

    Σ₁₁ = Σ₂₂ = m₀(4m₀+3)/((1+m₀)²(3m₀+2)) = 0.350218694…,      Σ₁₂ = −m₀(2m₀+1)/((1+m₀)²(3m₀+2)) = −0.148856754…,
    Σ₁₁ − Σ₁₂ = 2m₀/(1+m₀)²,   Σ₁₁ + Σ₁₂ = 2m₀/((1+m₀)(3m₀+2)),   ρ = Σ₁₂/Σ₁₁ = −(1+2m₀)/(3+4m₀) = −0.425039429596747….

*Proof.* (1) *Quadratic variation in natural time.* Let A^(ij)_k = Σ_(l<k) E[ΔLⁱ_(l+1) ΔLʲ_(l+1) | F_l]. By Proposition 2 (b),
A^(ij)_k = H_ij·T_k − β_ij·N(k) − γ_ij·Δ_k + R_k, where R is a martingale with bounded increments (the sum of the
differences between the three quantities and their conditional expectations). By Doob's inequality
max_(k≤n) |R_k| = O_P(√n), and by Lemma 4 max_(k≤n) |N(k)| = O_P(√n). So max_(k≤n) |A^(ij)_k − H_ij·T_k| = O_P(√n).
No ergodic theorem is needed: in natural time the quadratic variation is deterministic up to O(√n).
(2) *Martingale functional central limit theorem.* k(s) is a bounded stopping time, k(s+1) − k(s) ≤ 1 and
T_(k(s)) ∈ {s, s+1}. So Mⁿ_t = n^(−1/2)·L_(k(⌊nt⌋)) is a martingale with jumps ≤ B/√n, Mⁿ,ⁱMⁿ,ʲ − n^(−1)A^(ij)_(k(⌊nt⌋)) is
a martingale, and n^(−1)A^(ij)_(k(⌊nt⌋)) → H_ij·t in probability for each t, with jumps ≤ B²/n. By the martingale
functional central limit theorem (Ethier and Kurtz, Markov Processes, Theorem 7.1.4; cited from memory, not checked at
the book), Mⁿ converges to a Brownian motion with covariance H.
(3) *From L to the heights.* By (3.1) and Lemma 4, max_(k≤n) |(x_k, y_k)ᵀ − J⁻¹L_k| ≤ C·n^(0.4) + C with probability
1 − o(1). So the height process converges to J⁻¹·(Brownian motion with covariance H), with covariance
Σ = J⁻¹HJ⁻ᵀ. `stationary.py` computes Σ in Q(m₀) and verifies the four closed forms and ρ exactly. ∎

*Mixing form.* The limit in the martingale central limit theorem with a deterministic limiting quadratic variation
is mixing in the sense of Rényi (Hall and Heyde, Martingale Limit Theory, section 3.2; cited from memory), so
Theorem 3 holds also under P_π(· | A) for every event A of positive probability, in particular for every prescribed
level and top windows at time 0.

**Corollary 2 (the cone).** The linear map that turns Σ into the identity sends the quarter plane to a wedge of angle

    θ = arccos(−ρ) = arccos((1+2m₀)/(3+4m₀)) = 1.1317908814… = 64.846841 degrees.

For planar Brownian motion in a wedge of angle θ the survival exponent is π/(2θ) and the exponent of the excursion from
the apex back to the apex is 1 + π/θ. So

    α = 1 + π/arccos((1+2m₀)/(3+4m₀)) = 3.77577130646…

is, rigorously, the excursion exponent of the scaling limit of the bulk chain. This closes the gap named in 622
section 4 (iii): the covariance of 622 is the covariance of a proved limit theorem. What is not proved is that the
discrete return probability has this exponent (section 6).

## 6. What is missing for the exponent [not proved]

**Lemma L (open).** P_π(E_n | level 0, s₁ = −, t₁ = +) = Θ(n^(−α)) along even n, with α as in Corollary 2.

By (2.1), Lemma L is equivalent to the target M_m = Θ(Λ^m·m^(−α)).

**Why Theorem 2 of 227 does not apply.** Its hypotheses (A1) to (A4), (S1), (S2) are about a walk that is a sum of
i.i.d. cycles (W, Z, T) with an exponential moment (A1) and a killing rule that looks at one cycle and its starting
point (S1, S2). Here:
- Inside the excursion there is no time at which the future is independent of the past given a finite mode: every
  buried letter is a function of the past (section 1) and is used again when it comes back to the top.
- In the stationary bulk chain the only regeneration set found is R (simultaneous running minima, Corollary 1). Its
  gaps have no exponential moment, so (A1) fails, and the excursion never visits R after time 0.
- A hidden bijection with a finite-mode quarter-plane walk is not excluded. I looked for one in three ways and found
  none: (a) letters decided at the pop, fair or biased (section 1: fails exactly at m = 4); (b) coarser stack
  letters (blocks "− followed by a run of +"): the length of a buried run is again hidden information; (c) regeneration
  at running minima of one stack only: the other stack still carries the past.

**What is now available for a proof of Lemma L, and what a proof still needs.**
Available (all proved above):
1. Martingale coordinates: J⁻¹L is a two-dimensional martingale with bounded increments whose quadratic variation in
   natural time is Σ·T up to an explicit function of the state (Proposition 2 (b)): the usual corrector for the modes
   is explicit, and no mixing of the modes is needed.
2. The heights differ from J⁻¹L by the letter imbalances D^σ, D^τ only (3.1), and these are sums of the top letters
   (Lemma 3).
3. Under P_π the imbalances at depth j are of order √j with Gaussian tails at every time. By (2.1) and Theorem 2,
   P_π(E_n | tops) ≥ exp(−C√n), and it should be polynomial; an event of P_π-probability exp(−n^ε) is then negligible
   inside E_n as soon as a polynomial lower bound is known.

Still needed:
- **(L1) survival.** P_π(the initial tops are not popped during n points | tops) = Θ(n^(−π/(2θ))). In martingale
  coordinates this is the survival of a planar martingale with bounded increments and deterministic quadratic
  variation in a wedge whose boundary is displaced by the letter imbalance: at height y of the other stack the
  boundary x = 0 sits at distance of order √y from the boundary of the wedge. The method of Denisov and Wachtel (2015,
  Theorems 1 and 3) uses independent increments in the layer near the boundary; that part has to be redone.
- **(L2) gluing.** The return probability is a sum over the state s at half time of P(t → s)·P(s → t), and s is a pair
  of words, not a pair of heights. A local limit theorem for the heights is not enough: one needs that, given the
  heights, the law of the words produced by the first half and the law required by the second half are close to the
  uniform law in a sense strong enough to multiply. For walks with i.i.d. cycles this step is Lemma "exact
  decomposition" of 227 and the local limit theorem of Lawler and Limic (Theorem 2.3.5); here it is new.
- A weaker statement M_m = Λ^m·m^(−α+o(1)) needs the same two steps with losses n^(o(1)). I did not prove it. The
  published analog (Sheffield's hamburger-cheeseburger walk, where the increments are also read off a stack) has the
  exponent only in this weaker form.

## 7. The one-sided weighted family (weight t per point with g = 2) [family proved in 622; the rest data or conditional]

For the family of 622 section 5 the same structure holds at every t tried. `tfamily.py` (40 digits, t = 0.2, 0.6, 1, 2,
4) finds:
- the bulk critical chain has a product stationary law: the letters of σ are + with probability p_+q_−, those of τ with
  probability q_+p_− = 1 − p_+q_−, independent, except the active top; the level is not uniform for t ≠ 1 (P(level 1) =
  0.3678, 0.4383, 0.5, 0.5955, 0.6717);
- the first derivatives of the family are martingales, the heights have zero drift, and the covariance of the
  martingale increments per point equals minus the Hessian of log z (difference below 10⁻⁴⁰);
- the resulting correlation equals the formula of 622 (`family12_t.py`) to 15 digits.

So Theorems 1, 2, 3 extend to every t > 0 by the same proofs, once the finite check of Theorem 1 is done at that t.
Status: **data** for t ≠ 1 (the check was numerical, not symbolic in t). Given that check:
- Λ(t) = min over p of t(1+p)(1+μ)²/p, μ² = p(1+p)/(t(1+p) − p), is the exact growth rate of the weighted sum
  (upper bound proved in 622; lower bound by the proof of Theorem 2);
- the scaling limit of the bulk chain has correlation ρ(t), and α(t) = 1 + π/arccos(−ρ(t)) is its excursion exponent;
  that the weighted counts have this exponent is Lemma L at parameter t.

**Exact expansion at the Baxter point** (`tseries_exact.py`, rational power series; the coefficients below are exact
consequences of the closed formulas of the family and of Σ = −J⁻ᵀHJ⁻¹). With p = c·t the critical equation is
c(1+ct)(ctD + c)² = D³, D = 1 − c + ct, and

    c(t) = 1/2 + t/8 − t²/24 − 61t³/1152 − 31t⁴/3456 + …,        μ(t) = 1 − t²/12 − t³/36 + 29t⁴/1728 + …,
    Λ(t)/8 = 1 + t/4 + t²/16 + 0·t³ − 7t⁴/768 − 7t⁵/2304 + …,
    log(Λ(t)/8) = t/4 + t²/32 − t³/96 − 25t⁴/3072 − 13t⁵/23040 + …,
    ρ(t) = −1/2 + t/8 + 0·t² − 7t³/128 − 7t⁴/256 + 157t⁵/9216 + …,
    α(t) = 4 − (3√3/(4π))·t + ((18 + √3·π)/(32π²))·t² + O(t³).

The numbers 1/4, 1/32, −1/96, 1/8, 0, −7/128 and 0.0742222 of 622 (recognized there from 12 digits) are confirmed as
exact rationals, and the second coefficient of α gets the closed form (18 + √3π)/(32π²) = 0.074222193780746.
α′(0) = −3√3/(4π) follows from ρ′(0) = 1/8 alone: with θ = arccos(−ρ), θ(0) = π/3 and dθ/dρ = 1/sin θ = 2/√3,
dα/dt = −(π/θ²)·(2/√3)·(1/8) = −9/(4√3·π) = −3√3/(4π).
Status of the series for Λ: exact for the proved upper bound. Status of the series for ρ and α as statements about
the counts: conditional on the finite check at general t (data) and on Lemma L.

## 8. Algebraic facts [proved; `algebra.py`]

- 4x⁴ − 2x − 1 is irreducible over Q.
- Λ has minimal polynomial 16x⁴ − 112x³ − 520x² − 536x + 25 (irreducible): Λ is algebraic of degree 4. It is not an
  algebraic integer.
- c = cos θ = (1+2m₀)/(3+4m₀) has minimal polynomial 113x⁴ − 156x³ + 78x² − 16x + 1 (irreducible), and 2c has the
  primitive minimal polynomial 113x⁴ − 312x³ + 312x² − 128x + 16, which is not monic. So 2c is not an algebraic integer.
- **θ/π is irrational.** If θ = 2πk/n, then 2cos θ = ζ + ζ⁻¹ with ζ a root of unity, an algebraic integer. Contradiction.
- **θ/π and α are transcendental.** e^(iθ) = c + i·sqrt(1 − c²) is algebraic. If β = θ/π were algebraic, it would be
  algebraic and irrational, and (−1)^β = e^(iπβ) = e^(iθ) would be transcendental by the Gelfond–Schneider theorem.
  So β is transcendental, and so is α = 1 + 1/β.
- Consequence, if Lemma L holds: the generating function of M^[−1,2] is not D-finite (a D-finite series with integer
  coefficients and asymptotics C·Λ^m·m^(−α) has rational α; this is the argument of Bostan, Raschel and Salvy used in
  227). Without Lemma L nothing about D-finiteness is proved here.
