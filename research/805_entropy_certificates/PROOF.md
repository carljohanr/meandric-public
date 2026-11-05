# 805: entropy certificates (proofs)

2026-10-10. Status words as in `WAVE_620_COMMON.md`: **proved** (complete argument written here), **sketch**,
**conjecture**, **data**. Nothing here is refereed. Inputs from the project: Proposition T of
`637_rung2_matrix_product/PROOF.md` (tape machine; M_m = number of walks root → root of length 2m − 2), and the
ideas of `804_bulk_pressure/PROOF.md` section 4 (Theorem L, Lemmas R1, B1), which are re-proved below in the form
used, so that sections 1 and 2 depend only on Proposition T. M_m = number of meanders of order m, μ = lim M_m^(1/m).

Contents: 1 the stationary-flow theorem (Theorem F); 2 the exact certificate, μ ≥ 12.058577105; 3 Theorem L and the
recurrence lemma (question 1); 4 average against minimum (question 2); 5 explicit processes (question 3).

## 0. Conventions

Tape = word over {+, −} with a head on one letter. Moves and labels (637, 804): L, R (head one letter to the left,
right); splits S1: [+] → + [−] +, S2: [−] → + [−] −, S3: [−] → − [−] +; merge M: x [+] y → [x + y + 1], allowed
iff the head has both neighbors and (x, y) ≠ (+, +) (so − [+] − → [−], and the two mixed cases give [+]). Root:
+ [−] +. a, b = numbers of letters left and right of the head; a move changes (a, b) by L (−1, +1), R (+1, −1),
split (+1, +1), merge (−1, −1). On a finite tape the successors of a state under different labels are different
states, so walks and label sequences from a given tape correspond one to one. c_n = number of walks root → root
of length n; c_(2m−2) = M_m; closed walks concatenate, so c_n c_n′ ≤ c_(n+n′) and c_n ≤ μ^(n/2) for every n.
Class C: first and last letter +, sum of the letters 1; it contains the root and is closed under the moves.
Ω = {+, −}^ℤ is the set of two-sided infinite tapes with the head at place 0 (re-indexed after every move); all
moves act on Ω by the same local rules.

## 1. The stationary-flow theorem [proved]

### 1.1 Statement

**Data.** An integer r ≥ 2, n = 2r + 1.
- *Far tape.* Positive numbers b3(x, y, z), x, y, z ∈ {+, −}, with sum 1 and Σ_x b3(x, y, z) = Σ_x b3(y, z, x) =: b2(y, z)
  for all y, z. β is the stationary order-2 Markov measure on {+, −}^ℤ with these 3-gram probabilities. Put
  ρ(z | x, y) = b3(x, y, z)/b2(x, y) (the letter right of the pair x y) and λ(f | x, y) = b3(f, x, y)/b2(x, y) (the letter
  left of the pair x y).
- *Flows.* Numbers J(w, m) ≥ 0 for every window w ∈ {+, −}^n (the letters at the places −r, …, r) and every label m
  applicable at w (applicability depends on the three middle letters), with total 1. Put ν_W(w) = Σ_m J(w, m).

**The law ν and the chain.** ν is the probability measure on Ω under which the window has law ν_W and, given the
window w = w_1 … w_n, the letters to the right are produced one by one with the probabilities ρ(· | the two letters
before) starting from the pair w_(n−1) w_n, the letters to the left with λ(· | the two letters after) starting from
the pair w_1 w_2, independently on the two sides. From ω the chain makes the move m with probability
p(m | w) = J(w, m)/ν_W(w), w the window of ω.

**Conditions.**
- (S) *Stationarity.* For a word V = f2 f1 w g1 g2 (2r + 5 letters, head on the middle letter) put
  P(V) = ν_W(w)·λ(f1 | w_1, w_2)·λ(f2 | f1, w_1)·ρ(g1 | w_(n−1), w_n)·ρ(g2 | w_n, g1). For a label m applicable at w let
  W′(V, m) be the 2r + 3 letters around the head of the tape obtained from V by the move m (the head of the new
  tape always has at least r + 1 letters on each side). Required, for every W′ ∈ {+, −}^(2r+3) with inner window
  u = W′_2 … W′_(2r+2):

      Σ over (V, m) with W′(V, m) = W′ of  P(V)·p(m | w)   =   ν_W(u)·λ(W′_1 | u_1, u_2)·ρ(W′_(2r+3) | u_(n−1), u_n).

- (D) *Zero drift.* Σ_(w,m) J(w, m)·Δa(m) = 0 and Σ_(w,m) J(w, m)·Δb(m) = 0.
- (P) *Positivity.* J(w, L) > 0 for every window w.

**Theorem F.** If (S), (D), (P) hold, then

    μ ≥ exp(2h),        h = Σ_(w,m) J(w, m)·log(ν_W(w)/J(w, m))        (terms with J = 0 left out).

h is the entropy of the next move given the window, averaged over ν.

### 1.2 Proof

**Step 1: ν is invariant.** Let ω have law ν, choose m with the probabilities p(· | window of ω), and let ω′ be
the new tape. Let V be the 2r + 5 letters of ω around the head; under ν it has the law P(V), and given V the
letters right of V continue the chain from (g1, g2), those left of V from (f2, f1), independently; the label
depends on w only. Let W′ be the 2r + 3 letters of ω′ around its head; W′ = W′(V, m). In each case the letters of ω′
outside W′, given (V, m), are chain continuations from the two outer letters of W′ on each side:
- m = L: ω′_j = ω_(j−1), W′ = f2 f1 w. Right of W′ come g1, g2, …, the continuation from (w_(n−1), w_n), the last two
  letters of W′; left of W′ the continuation from (f2, f1), its first two letters. R is the mirror image.
- m a split: W′ = w_1 … w_r x′ − y′ w_(r+2) … w_n. Its first two letters are w_1 w_2 and its last two w_(n−1) w_n
  (r ≥ 2), and the far letters of ω are unchanged: f1, f2, … continue from (w_1, w_2), and g1, g2, … from (w_(n−1), w_n).
- m = M: W′ = f2 f1 w_1 … w_(r−1) z w_(r+3) … w_n g1 g2; the letters beyond continue from (f2, f1) and from (g1, g2).
So the law of ω′ is: W′ has the law given by the left side of (S), and the rest is continued by the chain from the
ends of W′. Under ν the 2r + 3 letters around the head have the law on the right side of (S) (window law times one
chain step on each side) and the rest is continued by the chain from their ends (the chain has order 2). Hence
(S) holds iff ω′ has law ν. ∎

**Step 2: the stationary chain is ergodic (Lemma G).** β is invariant and ergodic under the shift τ of {+, −}^ℤ
(an order-2 chain with positive transition probabilities). ν has the density ψ(window) = ν_W(w)/β(w) with respect to
β, and ψ > 0 everywhere by (P), so ν and β have the same null sets. Let g be bounded and measurable with
Σ_m p(m | ω) g(mω) = g(ω) for ν-almost every ω (a harmonic function). For the stationary chain (ω_0, ω_1),
E(g(ω_1) − g(ω_0))² = 2E g(ω_0)² − 2E g(ω_0)·(Pg)(ω_0) = 0 (invariance, then harmonicity). So g(mω) = g(ω) for
ν-almost every ω and every m with p(m | ω) > 0, in particular for m = L (by (P)), which is the shift:
g∘τ = g almost everywhere for ν, hence for β, hence g is constant. The same holds for the chain of pairs
(tape, next label): a bounded harmonic function G(ω, m) of that chain equals Ḡ(mω) with
Ḡ = Σ_m′ p(m′ | ·) G(·, m′) harmonic for the tape chain, hence constant. A stationary Markov chain whose bounded
harmonic functions are constant is ergodic (for a shift-invariant event E, the function ω ↦ P_ω(E) is harmonic,
and P_(ω_t)(E) → 1_E by martingale convergence, so 1_E is constant). ∎

**Step 3: two linear-time lemmas (804, Lemmas R1 and B1; proofs written out).**
A merge x [+] y → [x + y + 1] deletes the head letter + and one neighboring letter − and leaves the other neighbor
unchanged (− [+] − → [−]: delete the left −; − [+] + → [+]: delete the −; + [+] − → [+]: delete the −). A split read
backward, x′ [−] y′ → [x′ + y′ − 1] with (x′, y′) ≠ (−, −), deletes the head letter − and one neighboring letter + and
leaves the other neighbor unchanged. In both cases the head has both neighbors, so the two deleted letters are
not the first or the last letter of the tape when the deleted neighbor is an inner letter.

*Lemma R (return).* From a tape of class C with N letters the root is reached in at most 3N moves through tapes
with at most N letters. *Lemma B (build).* Every tape of class C with N letters is reached from the root in at
most 3N moves through tapes with at most N letters.

*Proof.* Number the places 0, …, N − 1; the inner places are 1, …, N − 2; the inner word has sum −1. Sweep: put
the head on place 1 (at most N moves) and repeat, with the head on place h, keeping the invariant "the inner
letters at places < h are equal":
 (i) if h and h − 1 are inner places carrying different letters, delete this pair: for Lemma R by a merge with the
 head on the + of the pair (zero or one head move first), for Lemma B by a backward split with the head on the −
 of the pair (zero or one head move first). The deleted letters are inner, so the word stays in class C; the head
 ends on place h − 1 or h − 2 and the invariant holds;
 (ii) otherwise move right; the invariant holds because the letters at h − 1 and h are equal or h ≤ 1.
If the head reaches the last place, all inner letters are equal with sum −1, so the word is + − +. There are
(N − 3)/2 deletions, each with at most 2 moves and a decrease of h by at most 2, so at most N − 1 moves right; with
the first positioning and one last move onto the −: at most 3N moves, and the length never grows. For Lemma B
read the walk backward: head moves are reversible and a backward split read forward is a split, which is always
allowed. ∎ (804 checked both procedures on all 33,535 headed tapes of class C with at most 7 pairs.)

**Step 4: counting (804 Theorem L in the Markov case).** Run the stationary chain: ω_0 with law ν, labels
m_1, m_2, …, A_i and B_i the total change of a and b after i moves, I_n = −Σ_(i ≤ n) log p(m_i | ω_(i−1)).
By Step 2 and the ergodic theorem, A_n/n → Σ J·Δa = 0, B_n/n → 0 (condition (D)) and I_n/n → h almost surely
(−log p(m_1 | ω_0) is non-negative with mean h ≤ log 4). Fix ε > 0. For large n the event
E = {abs(A_i) + abs(B_i) ≤ εn for all i ≤ n, and I_n ≥ n(h − ε)} has probability at least 1/2, so there is an ω with
P(E | ω_0 = ω) ≥ 1/2. Given ω_0 = ω a label sequence has probability exp(−I_n) ≤ exp(−n(h − ε)) on E, so E contains
at least (1/2)·exp(n(h − ε)) label sequences, each applicable from ω (a label of positive probability is applicable).
Let ℓ = ⌈εn⌉ + 2 and let T_0 be the 2ℓ + 1 letters of ω around the head. Along each of these sequences the head keeps
at least two letters of T_0 on each side, and moves read and change only the head letter and its two neighbors,
so each sequence is a walk of the tape machine from T_0. Put + c^i in front and c^(i′) + behind, with one sign c
and i + i′ = abs(1 + sum of T_0), so that the result T_1 is of class C; it has at most 4ℓ + 5 letters and the same
walks. Go from the root to T_1 (Lemma B, one fixed walk of at most 3(4ℓ + 5) moves), follow one of the walks, and
return to the root from its end tape, which is of class C with at most 4ℓ + 5 + εn letters (Lemma R). Different
label sequences give different closed walks from the root (they differ inside the middle part), with lengths in
[n, n + n′], n′ ≤ 30εn + 110. So for some N in this range c_N ≥ exp(n(h − ε))/(2(n′ + 1)), and c_N ≤ μ^(N/2) gives
n(h − ε) − log(2n′ + 2) ≤ ((n + n′)/2)·log μ. Divide by n, let n → ∞, then ε → 0: h ≤ (1/2) log μ. ∎

Remarks. 1. Theorem F is Theorem L of 804 for a Markov chain whose stationary law is known in closed form;
what is new is the class of laws (window law times a Markov far tape), the reduction of stationarity to the
finite linear system (S), and Lemma G, which supplies the ergodicity that 804 left open for product laws.
2. For a *fixed* far-tape law, (S), (D) and the normalization are linear in J and h is a concave function of J.
So the best process of this shape is the solution of a finite concave maximization over a polytope, and every
feasible point is a bound. This is what makes the format a certificate format.
3. On a strip of K levels the same statement and proof hold with a nearest-neighbor Markov path as far tape
(`flow.py`). For the strip of three the maximum is 10.4245356612 at radius 1, 2 and 3 (data, 10 digits), the exact
growth rate of that band (625, 634): the program finds the product law of 634 and nothing above it.

## 2. The exact certificate: μ ≥ 12.058577105 [certified (exact computation), given Theorem F]

**Certificate** `cert_r5.pkl` (879 KB; sha256 in `SHA256SUMS`): r = 5; far tape = the order-2 chain with
P(+ | − −) = s = 23/32, P(+ | − +) = t = 7/16, P(+ | + −) = 1 − t, P(+ | + +) = 1 − s (pair law 7/37 for equal letters,
23/74 for different letters); 7,936 flows J(w, m) = (integer)/(common denominator of 824 bits) for the 2,048
windows of 11 letters.

**Verifier** `verify_cert.py` (written from the definitions in section 1.1 and the machine of 637/641; it shares
no code with the generator). It checks, in exact integer arithmetic: the far-tape law is a stationary 3-gram law
with positive entries; the flows are non-negative integers, one for every applicable (window, label) and no
others, every L-flow is positive, the total is the denominator; (S) by brute force: for each of the 2,048 × 16
words V of 15 letters it makes every applicable move on the word with the machine, reads the 13 letters around
the new head and accumulates the integer weight, then compares with the right side for all 8,192 windows W′;
(D). Then it evaluates h by interval arithmetic (mpmath.iv, 200 bits). It first tests its own machine: the
closed walks from the root number 1, 2, 8, 42, 262, 1828, 13820. Output (`log_verify_r5.txt`, 0.5 s):

    (V2) stationarity: all 8192 equations hold exactly        (V3) zero drift of both side lengths: exact
    h ≥ 1.24488809985726960773743460491      exp(h) ≥ 3.47254619917969857045714287994
    exp(2h) ≥ 12.0585771054373707766456619765         exp(2h) − (3471821951/10⁹)² ≥ 0.00502944599192
    ACCEPTED: mu >= 12058577105/10^9

**Corollary 2.1 (given Theorem F and this computation).** μ ≥ 12.058577105 (√μ ≥ 3.4725461991). This is above the
record of 603, μ ≥ 12.053547659445446401, by 0.00503.

Teeth (`teeth.py`, `log_teeth_r5.txt`): seven damaged copies (one flow + 1; L + 1 and R − 1 at one window, which
keeps the window law and the total; S2 and S3 exchanged at one window; s changed; one L-flow set to 0;
denominator + 1; one flow missing) are all rejected. A second certificate `cert_r3.pkl` (r = 3, 496 flows, 15 KB,
denominator of 175 bits) is accepted with μ ≥ 11.843448806.

How the certificate was made (`exact_cert.py`; not part of the proof): float maximum of the concave program
(`flow2.py`, Newton steps on the sparse optimality system), made mirror symmetric; exact integer rows on the
mirror orbits (3,992 unknowns, rank 2,791); 1,201 free flows rounded to multiples of 2⁻⁵⁶; the other 2,791 from
the square system, solved exactly by refinement with integer residuals and one rational reconstruction
(denominator of 768 bits); 50 s in all.

What the bound rests on: Proposition T of 637 (reviewed in 677); Theorem F (section 1, written here, not
refereed; Step 4 is Theorem L of 804, also not refereed); the verifier and its run. One author wrote generator
and verifier. Independent evidence that the method does not overcount: on the strip of three the program returns
the exact growth rate 10.4245356612 at radius 1, 2, 3 and not more; on the strip of four it gives 11.2692
(radius 2) and 11.3666 (radius 3), below the growth rate 11.404 of rung 2; every value for all meanders is below
the estimate 12.2629 of μ and the proved upper bound 12.4147.

**Float values (data; `DATA.txt`).** exp(2h) for all meanders, far tape independent fair letters / best order-1
chain / order-2 chain: radius 1: 10.9976; radius 2: 11.4946 / 11.5567; radius 3: 11.7345 / 11.8310 / 11.8434;
radius 4: 11.8409 / 11.9464 / 11.9875; radius 5: 11.9151 / 12.0173 / 12.0590; radius 6, order 1: 12.0589.
Number of flows 3.875·2^(2r+1): 31, 124, 496, 1,984, 7,936, 31,744. For comparison, the indexed bracket table of
603 has 3.4210 = √11.703 at radius 3 with 0.38 million entries (J = 80) and needs 1.26 billion entries for
3.4718; here 496 flows give 3.4414 and 7,936 flows give 3.4725.

## 3. Theorem L and the recurrence lemma (question 1)

Setting for this section: S = the tapes of class C with a head (a countable set), A the 0/1 matrix of the moves,
o the root. u ≥ 0 is a weight with support D ∋ o and (Au)(s) ≥ R·u(s) for s ∈ D, R > 0 (a *sub-solution*; A u is
the sum over successors, and successors outside D count 0). R′(s) = (Au)(s)/u(s) ≥ R and
P(s, s′) = u(s′)/(Au)(s) is the transformed chain on D. All statements hold verbatim for the bracket-word machine.

**Lemma 3.1 (spectral form of the recurrence lemma) [proved].** For every n,
√μ ≥ R·P^n(o, o)^(1/n). Hence √μ ≥ R·ρ(P), ρ(P) = limsup P^n(o, o)^(1/n), and μ ≥ R² as soon as the return
probabilities of the transformed chain do not decay exponentially. Recurrence (Σ_n P^n(o, o) = ∞) is sufficient
and not necessary.

*Proof.* A walk s_0 → … → s_n in D has P-probability (u(s_n)/u(s_0))·Π_(i<n) 1/R′(s_i). For closed walks from o
this is at most R^(−n), so P^n(o, o) ≤ R^(−n)·c_n ≤ R^(−n) μ^(n/2). ∎

**Lemma 3.2 (increments of a stationary sequence) [proved; standard].** Let (X_t) be stationary and F a
measurable real function such that the negative part of D = F(X_1) − F(X_0) is integrable. Then D is integrable
and E[D | invariant σ-field] = 0; in particular E[D] = 0.

*Proof.* Birkhoff's theorem for functions with integrable negative part: (1/n)(F(X_n) − F(X_0)) → E[D | 𝓘] almost
surely, with values in (−∞, +∞]. The left side tends to 0 in probability because F(X_n) has the law of F(X_0). ∎

Consequences. (a) A stationary process on *finite* tapes has zero mean change of a and of b: zero drift is
automatic there and is a genuine hypothesis only for infinite tapes. (b) If P is positive recurrent with
stationary law π, then E_π[log R′] is finite and equals the entropy rate h(P) = E_π[−Σ_s′ P(s, s′) log P(s, s′)]:
apply the lemma to F = log u, D = log P(X_0, X_1) + log R′(X_0) ≥ log P(X_0, X_1) + log R, whose negative part is
integrable because −log P(X_0, X_1) ≥ 0 has mean at most log 4.

**Proposition 3.3 (average form; Theorem L on finite tapes with an elementary proof) [proved].** Let p be any
Markov kernel on a set D ∋ o of tapes of class C, supported on the moves, irreducible and positive recurrent
with stationary law π. Then log √μ ≥ h(p) = E_π[−Σ p log p]. For p = the transformed chain of a sub-solution u:

    log √μ  ≥  E_π[log R′]  ≥  log R,         with equality on the right iff R′ = R on the support of π.

*Proof.* With I_n = −Σ_(i<n) log p(X_i, X_(i+1)) ≥ 0 and E_o the chain started at o,
c_n ≥ Σ over closed walks in D of 1 = E_o[exp(I_n); X_n = o] ≥ P^n(o, o)·exp(E_o[I_n | X_n = o]) (Jensen).
For K > 0 let I_n^K be the same sum with each term cut at K. By the ergodic theorem for the chain of pairs
(X_i, X_(i+1)) started at o, I_n^K/n → h_K = E_π[Σ p·min(−log p, K)] almost surely and boundedly, so
E_o[I_n^K/n; X_n = o] − h_K·P^n(o, o) → 0; along the multiples of the period, P^n(o, o) tends to a positive limit.
Hence liminf (1/n) log c_n ≥ h_K, and h_K ↑ h(p). ∎

This is the easy half of Gurevich's variational principle for countable Markov shifts, in the Markov case. It
uses neither Lemma R nor Lemma B nor the Shannon–McMillan theorem.

**How the two tools relate (exact hypotheses).**
1. *Theorem L ⇒ recurrence lemma, with a better constant, when the transformed chain is positive recurrent:*
   it is then a stationary ergodic process on finite tapes, zero drift is automatic (3.2 (a)), its entropy rate
   is E_π[log R′] (3.2 (b)), and Theorem L (or 3.3) gives μ ≥ exp(2E_π[log R′]) ≥ R².
2. *Null recurrent or transient chains:* Theorem L does not apply (there is no stationary probability law);
   Lemma 3.1 does, and it needs only ρ(P) = 1. Through finite sections Theorem L still implies it: for finite
   F ∋ o, the Parry measure of the class of o in A restricted to F is an ergodic process with entropy
   log ρ(A_F) ≥ log(R·ρ(P_F)) (3.1 applied inside F), and sup_F ρ(P_F) = ρ(P) because every P^n(o, o) is a finite sum.
3. *Recurrence lemma ⇒ Theorem L only for chains of gradient form* P(s, s′) ∝ u(s′), on finite tapes, and then
   with the minimum of R′ in place of the average. For the average one needs a corrector (section 4) and
   recurrence of the *corrected* chain, which is automatic only for finite D. For processes with an infinite far
   tape (the case of section 1) there is no reference state and no weight u; Theorem L is the only tool.
4. So each proves the other exactly in the intersection "positive recurrent transformed chain of a
   sub-solution on finite tapes", where 3.1 and 3.3 give two independent short proofs.

**Proposition 3.4 (the recurrence lemma from the finite-support theorem, under Foster's criterion) [proved].**
Let u be a sub-solution with constant R on a possibly infinite support D, and let L: D → [0, ∞) have finite
sublevel sets and satisfy, for a finite set F ⊂ D,

    Σ_s′ u(s′)·(L(s′) − L(s)) ≤ 0   for s ∈ D \ F        (L is superharmonic for the transformed chain off F).

Then for every small c > 0 the weight w_c = u·(1 − cL)⁺ has finite support {L < 1/c} and is a sub-solution with
constant R_c = R·(1 − c·δ/(1 − c·max_F L)), δ = max over F of (PL − L)⁺. Hence μ ≥ R² by the finite-support
theorem (601, 603), with no Markov-chain theory.

*Proof.* For s with L(s) < 1/c: (A w_c)(s) ≥ Σ_s′ u(s′)(1 − cL(s′)) = (Au)(s)·(1 − c(PL)(s)) ≥ R′(s)·u(s)·(1 − cL(s) − c(PL − L)⁺(s)).
Off F this is R′(s) w_c(s). On F, u(s) ≤ w_c(s)/(1 − c·max_F L). Let c → 0. ∎

In the profile format of 801/802 both u and u·L are of the form (exact global factors)·θ^(i+m)·(polynomial in
(i, m) with window-table coefficients), so the two conditions (sub-solution; drift of L) are row inequalities of
the same kind, and the soundness proof is the three lines above.

**Bracket-word version of Theorem L [proved, with the two lemmas below as sketch checked by enumeration].**
State (s, i), labels OO, OC, CO, CC; (i, m) = symbols below and above the marker change by (+1, +1), (−1, +1),
(+1, −1), (−1, −1). What has to be checked, and how:
(1) *Entropy.* A label and the state before determine the state after, but not conversely (the far flip of a CC
is not determined by the result). Reversibility was used in 804 only to identify σ-fields; the proof goes through
with h defined as the entropy of the next label given the whole past trajectory, and for a Markov kernel on
states h = E_π[−Σ p log p] in any case.
(2) *Locality.* A walk of n labels from a two-sided infinite word depends on the far part only through far
flips. Cut a window of ℓ > εn + 2 symbols on each side and complete it to a balanced word in any way: a symbol next
to the marker has its partner inside the window in the infinite word iff it has it inside the window in the
completed word, and then at the same place; so the two evolutions agree inside the window as long as the marker
keeps one window symbol on each side. (In the infinite word a flip with no partner does nothing.)
(3) *Linear build and return* for balanced words with 2k symbols. Return: sweep the marker to the right and
delete every ")(" met (a CC there changes nothing else); what is left is (^j )^j; put the marker before the last
"(" and apply CC j − 1 times (each gives (^(j−1) )^(j−1) with the marker again before the last "("). Build: a
balanced word is made from "()" by insertions of "()" (labels OO), and the leftmost "()" can always be removed
last-in-first-out with total marker travel at most 3k. Both are O(k) moves [sketch; not machine-checked here].
(4) *Drift variables*: (i, m). With (1) to (4) the proof of Step 4 in section 1.2 applies word for word. For
processes on *finite* words Proposition 3.3 already gives the statement without (2) and (3).

## 4. Average against minimum (question 2)

**Format.** A *format* is: a finite or countable set of classes c(s) of states; rows ρ(s) such that the row fixes
the class c(ρ) and, for every move m, the class c′(ρ, m) of the successor and a factor g(ρ, m) > 0 (the exact
ratio of the global factors; a far flip with two possible peak changes is two rows). Weights of the format are
u(s) = G(s)·U(c(s)); then R′(s) = R′_U(ρ) = Σ_m g(ρ, m) U(c′(ρ, m))/U(c(ρ)) and the transformed chain makes the move m
with probability p^U_ρ(m) = g(ρ, m) U(c′(ρ, m))/(R′_U(ρ) U(c(ρ))). A *corrector* is a function φ of the class. Put

    point(U) = inf_ρ log R′_U(ρ),        avg(U, φ) = inf_ρ [ log R′_U(ρ) + Σ_m p^U_ρ(m)·φ(c′(ρ, m)) − φ(c(ρ)) ].

**Theorem D [proved].** (a) For every stationary law π of the transformed chain of U (and every law on states
whose class marginal is balanced), E_π[log R′_U] ≥ avg(U, φ): this is the bound on the average that linear
programming duality gives from the stationarity constraints with φ as dual variable, and for finitely many rows
sup_φ avg(U, φ) = min{Σ_ρ q(ρ) log R′_U(ρ) : q a probability vector on rows, balanced on classes}.
(b) For every U and φ:   point(U·e^φ) ≥ avg(U, φ), with equality iff φ∘c′ is constant on the moves of the binding
rows. (c) Consequently sup_(U,φ) avg(U, φ) = sup_U point(U): **a bound on the average obtained from stationarity
constraints with correctors in the format is never better than the pointwise sub-solution constant of the
corrected weight U·e^φ, which lies in the same format.** The claim of the brief is true.

*Proof.* (a) E_π[φ(c(X_1)) − φ(c(X_0))] = 0 (Lemma 3.2; trivial for bounded φ), so
E_π[log R′_U] = E_π[log R′_U + Pφ − φ] ≥ the infimum over rows. The second statement is the duality theorem of
linear programming (the primal is feasible and bounded). (b) R′_(Ue^φ)(ρ) = Σ_m g U(c′) e^(φ(c′))/(U(c) e^(φ(c)))
= R′_U(ρ)·Σ_m p^U_ρ(m) e^(φ(c′(ρ,m)) − φ(c(ρ))) ≥ R′_U(ρ)·exp(Σ_m p^U_ρ(m) φ(c′(ρ, m)) − φ(c(ρ))) by Jensen. (c) "≤" is (b), "≥" is
φ = 0. ∎

Remarks.
1. *Side conditions differ.* The averaged bound is a bound on μ when the chain of U is positive recurrent
(Proposition 3.3); the pointwise bound for U·e^φ needs ρ(P) = 1 for the chain of U·e^φ (Lemma 3.1), which is a
different chain. So (c) compares constants, and each constant needs its own recurrence statement, except with a
finite support, where both are automatic.
2. *Dictionary: every stationarity constraint is a factor.* For any function Φ of the state whose increment is
fixed by the row, E_π[ΔΦ] = 0 (Lemma 3.2), and the multiplier of this constraint is the factor e^(κΦ) of the
weight. Balance of window classes ↔ the window table; zero drift of the side lengths ↔ θ^(i+m) (804: length tilt);
zero flux of k-grams ↔ the k-gram factor of 641; the peak count ↔ p^(peaks); stationarity of the (i, m)-marginal ↔
the index (i, m) of 603, of which a polynomial profile is a subfamily. Nothing of this kind is gained by averaging.
3. *What the pointwise bound does not see.* Let 𝓜 be the set of balanced row laws q (the minimum in (a)). The
true stationary law π of the chain started at the root is one point of 𝓜; the others are stationary laws of a
game in which an adversary supplies the hidden letters at random, subject only to balance. The minimum over 𝓜 is
what a table certifies; E_π[log R′_U] can be larger. To certify more one needs constraints on π that are not of
the form E[ΔΦ] = 0: (i) *the far tape was written by the process itself*, which singles out π among the invariant
laws of the bulk chain (the sea of 800 is another invariant law of the same row dynamics, and it is what caps
the tape tables); with tables this can only be approached by indices that grow (radius, (i, m)), since each
radius adds its own balance constraints and leaves the outermost letter to the adversary; (ii) *knowing π
exactly*. Section 1 is (ii): the far tape is not left to an adversary but prescribed (a Markov chain), and (S) is
the condition that the process reproduces it. Then the average is exact, and no minimum is taken at all.
4. *What averaging adds in section 1.* In Theorem F the bound is an average, h = E_ν[log(1/p)]; no pointwise
inequality holds, and the kernel is not of gradient form. Its dual (the multipliers of (S), (D)) is an *upper*
bound certificate for the family of processes with that far tape, of geometric-mean type:
h_max = min over φ on windows of radius r + 1 and tilts of max_w log Σ_m t^(Δ(m)) exp(E_fresh[φ(W′)] − E_outer[φ(W ⊃ w)]).

## 5. Explicit processes with computable entropy (question 3)

**(a) Regenerative constructions.** Independent blocks between returns to the root: the bound is 1/x_N² with
Σ_(n ≤ N) f_n x_N^n = 1, f_n = number of first-return walks of length n, 1 − 1/Σ c_n x^n = Σ f_n x^n. Because
Σ c_n μ^(−n/2) < ∞ (the exponent of M_m is above 2 for every candidate value), Σ f_n μ^(−n/2) < 1, so x_N stays
away from μ^(−1/2) until the truncated tail is made up by a larger x: log(x_N √μ) ≈ (α log N)/N [sketch]. With
N ≤ 60 (the known M_m) this is of order 20 percent in √μ: far below 12 (not computed here; the value is of the
size of the old renewal bounds, ledger I36, J36). Regeneration at a bulk reference tape (for example the tape
of all −) is not easier: that tape reaches the padded root in five moves (S2, R, S2, L, L) and comes back in four,
so its closed bulk walks contain a copy of every walk root → root, and counting its return blocks is at least as
hard as counting meanders. Markov-regenerative
versions are finite graphs of tapes (Perron roots P_j of 800: 3.2912 at 9 pairs) or the window tables. So
(a) cannot approach 12.0535 with enumerable block sets.

**(b) The strip of three with excursions.** *Loop insertion [sketch; the exact form is Theorem 6 of 634].* Let ν be an ergodic zero-drift process
with entropy h_0 and let E be a set of tapes with ν(E) = e > 0 at which a closed walk γ of length ℓ (same tape and
same side lengths at the end) is applicable whose first label has ν-probability 0 at tapes of E. Insert γ with
probability ε at each visit to E. The new process is stationary after the time change, ergodic, has zero drift,
and its entropy per move is at least (h_0 + e·(H_2(ε) − ε·h_E))/(1 + εeℓ), H_2 the binary entropy and h_E ≤ log 4 the
mean entropy of the base move on E; this is h_0 + e·ε log(1/ε) + O(ε) > h_0 for small ε. So the product law of the
strip is never optimal among all meanders, the optimal ε is about exp(−ℓ h_0), and the gain is of order
e·exp(−ℓ h_0) per loop type. Summing all loops that leave the strip at an extreme level and come back to the same
tape is exactly the sealed-gadget substitution of 634 (Theorem 6, AU4, AU5): growth 1/x with x·Λ(S⁻(x), S⁺(x)) = 1,
values 10.7060 (S1), 10.9615 (S2), 11.0237 (Z2), 11.0331 (B2), ≥ 11.0462 (W2), and "towers add little". These are
instances of Theorem L (the base law is ergodic by Lemma G: the head moves have positive probability and the
law has a window density with respect to a shift-ergodic path measure; this closes the item "ergodicity of the
634 product law" of 804, section 6.6, up to rewriting 634 Theorem 5 (b) in tape coordinates [sketch]). Ceiling:
about 11.05 for loops at the extreme levels, 11.39 for the widest radius-1 product eigenfunctions (AU7, data).
(b) does not approach 12.0535.

**(c) Markov far tape, window-driven moves.** This is section 1. Stationarity is the linear system (S) once the
far-tape law is fixed; the far-tape law enters through at most 2^k numbers for order k. Differences to 634/637:
no right eigenfunction is needed (the kernel is free, not a Doob transform), so the no-go theorem of 637 for
matrix-product eigenfunctions does not apply, and the window may have any radius. Values: section 2. It passes
12.0535 at radius 5 with an order-2 far tape (certified) and at radius 6 with an order-1 far tape (float).
Gains seen (data): per radius step, in exp(h), 0.0741, 0.0352, 0.0155, 0.0108 (independent letters, radius 1 to 5)
and 0.0401, 0.0167, 0.0102, 0.0060 (order 1, radius 2 to 6); per order of the far tape at radius 4: 0.0153 (order
1), 0.0059 (order 2); at radius 5: 0.0148, 0.0059. Extrapolation [conjecture]: order 2 at radius 6 about 3.479
(μ ≥ 12.10), radius 7 about 3.483 (12.13). Whether the supremum over radius and order is √μ is open; it is at
most √μ by Theorem F and it would follow from an approximation of the Parry measures of 804 Theorem E by
processes with Markov far tape, which I did not attempt.
