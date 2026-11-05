# 1252 PROOF: the theory behind mu >= exp(2h) on a strip with a finite-memory reference law (my own proofs)

Agent EY, 2026-10-12. Status words as in `WAVE_1200_COMMON.md`. Everything here is **proved** (complete argument, written below) unless marked
otherwise. These proofs were written after reading the project's arguments and are not a copy of them; where they agree with 805 / 808 / 1008 /
1211 I say so. Numerical confirmations of each statement are in `FINAL.md` section 3 (scripts and logs in this folder).

Inputs taken as given: Proposition T of `637_rung2_matrix_product/PROOF.md` (reviewed in `677_tape_bijection_review`; I also re-counted it, see
FINAL): M_m = number of closed walks of length 2m - 2 from the root `+[-]+` of the tape machine; hence c_N := number of closed walks of length N
satisfies c_(2m-2) = M_m, c_N = 0 for N odd, c_N c_N' <= c_(N+N'), and by Fekete c_N <= mu^(N/2) for all N.

## 0. Setting

*Tape machine* (805 section 0). A tape is a finite word over {+,-} (a path of +-1 steps from a start height) with a head on one letter. Moves:
L, R (shift the head), splits S1 [+] -> + [-] +, S2 [-] -> + [-] -, S3 [-] -> - [-] +, and the merge x [+] y -> [x+y+1] (head + with both
neighbors, forbidden iff x = y = +). Only the local rules matter: a move reads and changes the head letter and its neighbors.

*Strip.* Fix K >= 2. A *level-labelled* tape has its start height; legal means all vertex heights lie in {0..K-1}. A split is legal iff the new path
stays in the strip (S1 always; S2 needs h+1 <= K-1; S3 needs h-2 >= 0 for a letter going from height h to h-1); merges and shifts keep the path.
K = infinity (no strip) is allowed everywhere below.

*Reference law beta of memory m.* A stationary Markov chain on the strip whose state is the last m edges (m+1 consecutive heights); assume
(H1) every admissible path of m+2 heights has positive probability, (H2) the state graph is irreducible, (H3) the transition kernel is stationary for
the block law. (The project's law is the quintet law t' with m = 3; I verified H1-H3 exactly for K = 12, kappa = 4/5, eta = 3/16, zeta = 1/8,
`t6_law_hypotheses.py`.) Reversibility is not used. beta(w) is the cylinder probability of a level-labelled word w (m+1 heights at least).
For a window W containing at least m letters at each end, beta(outside | W) is the product of the right chain continuation from the last m letters and
the left chain continuation from the first m letters (independent), and it depends only on those end blocks (and their heights).

*Windows, flows.* r >= 1, n = 2r+1. A window w is a level-labelled word of n letters with the head on the middle letter. A flow is a real function
J(w, l) on legal pairs (window, label l applicable at w). Put nu_J(w) = sum_l J(w,l). Rows: for every legal level-labelled word W' of n+2 letters with
inner window u (drop one letter each side),

    (S_r)   sum over (V,l) with W'(V,l) = W'  of  J(w,l) * beta(V)/beta(w)   [splits: J(w,l) * 1]   =   sum_l' J(u,l') * beta(W')/beta(u).

Here V runs over the words of n+4 letters whose middle n letters are w, and W'(V,l) is the (n+2)-letter word around the new head after l: for L
W' = V[0:n+2], for R W' = V[2:n+4]; for a split W' is w with the head letter replaced by three letters (no V needed); for the merge W' is V with the
three letters x,+,y replaced by one (n+4 -> n+2 letters). (This is exactly (S) of 805 written in terms of cylinder probabilities.)
Drifts (D): sum J * Delta a = sum J * Delta b = 0 with Delta(L) = (-1,+1), Delta(R) = (+1,-1), splits (+1,+1), merge (-1,-1). Mass: sum J = 1 when J >= 0.

*The measure of a flow.* Omega = bi-infinite level-labelled tapes (all admissible height paths in the strip) with a marked head; beta is a measure on
Omega (the stationary chain). For a flow J set mu_J(d omega, l) = J(w(omega), l)/beta(w(omega)) beta(d omega), w(omega) the radius-r window. mu_J is a
signed measure when J is signed; nu_J = sum_l mu_J(., l). T_l denotes the map of Omega given by the move l (defined where l is applicable).

## 1. Proposition A (finite rows give invariance of the whole tape law; r >= m) [proved]

**Statement.** Let r >= m and let J be any real flow satisfying (S_r). Then sum_l (T_l)_# mu_J(., l) = nu_J as signed measures on Omega.
Conversely invariance implies (S_r). In particular a flow satisfying (S_r) for r >= m also satisfies (S_R) for every R >= r after the lift
J_R(W,l) = J(center_r(W), l) beta(W)/beta(center_r(W)) (the lift is the same measure mu_J read at radius R).

**Proof.** Let C(W') denote the cylinder of tapes whose radius-(r+1) window is W'. Both sides are measures of the form
sum_W' mass(W') beta(. | W') over this partition once we show it for the left side; the right side is of that form because nu_J has density
sum_l' J(u,l')/beta(u) depending on the inner window u of W' only, so nu_J(C(W')) = RHS(W') of (S_r) and nu_J( . | C(W')) = beta(. | W').
Left side. Fix (V,l) with source window w and let Sset be the set of letters of omega that is read-or-changed together with what is kept in the target
window: Sset = W' for L (the letters f2 f1 w), Sset = V' = w g1 g2 for R, Sset = w for a split, Sset = V for a merge. On C(Sset) the measure
mu_J(., l) restricted is J(w,l)/beta(w) beta|_C(Sset); T_l changes only letters inside Sset, and maps the cylinder of Sset onto the cylinder of W'
with the letters outside Sset unchanged and in the same order. Condition E: the first m and last m letters of W' (with heights) are those of Sset.
Check of E: L, R: W' and Sset share their outer end blocks (for L, W' = Sset); split: the replaced letter is at index r+1 of w and indices m+1.. are
untouched, so the first m letters of W' are w_1..w_m (needs r >= m) and symmetrically at the right end; merge: the merged triple sits at indices
r+3..r+5 of V (1-based, V has 2r+5 letters), so the first m letters of V survive as soon as r + 2 > m. So E holds for r >= m (the split is the binding case).
By E and the Markov property of order m, beta(. | Sset) = beta(. | ends of Sset) = beta(. | ends of W') = beta(. | W') on the outside. Hence
(T_l)_# [mu_J(., l)|_C(Sset)] = J(w,l) beta(Sset)/beta(w) * beta(. | W'), a multiple of the same kernel. Summing over all (V,l) leading to W' gives
mass(W') = LHS(W') of (S_r) (beta(Sset)/beta(w) = beta(V)/beta(w) for merge, beta(W')/beta(w) for L, R, 1 for splits; for the merge several V can
lead to the same W', they all add). So the pushed measure equals sum_W' LHS(W') beta(.|W') and (S_r) says LHS = RHS. The converse is the same
computation read backwards. The lift statement: mu_J is one measure, its radius-R description is the lifted flow. QED

**Sharpness (data).** The proof needs r >= m and the numerics below show the statement fails for r = m - 1: for K = 3, 4, 5 and memory m = 2, 3 the
null space of the rows at radius r - 1 is strictly larger than the null space of the stacked system (rows at r-1 and lifted rows at r+1, r+2), and
equal for r >= m (FINAL section 3.2). So the project's hypothesis r >= m + 1 (1019) is sufficient and not sharp; nothing is lost at the used radii 6 to 10.
The proposition holds for signed J, which is what the signed-parent construction needs.

## 2. Theorem F' (strip, memory m, no strip needed for the completion) [proved]

**Statement.** Let beta satisfy H1-H3 with memory m on the strip of K levels (K = infinity allowed), r >= m, and let J >= 0 be a flow with (S_r), (D),
total mass 1 and J(w, L) > 0 for every legal window w (P). Let h = sum J log(nu_J/J) (terms with J = 0 omitted). Then mu >= exp(2h).

**Proof.** The proof is that of 805 / 808 with three points made explicit.
(1) *Invariance* is Proposition A: nu_J is a probability measure (J >= 0, mass 1), p(l | omega) = J(w,l)/nu_J(w) defines a Markov chain on Omega with
invariant law nu_J.
(2) *Ergodicity (Lemma G, strip version).* beta is the stationary law of an irreducible finite chain (H2), so the shift tau is ergodic for beta (period 2
does not matter: ergodic, not mixing). By (P) the density psi(w) = nu_J(w)/beta(w) is > 0 on every legal window, so nu_J and beta are equivalent. For a
bounded harmonic g (sum_l p(l|omega) g(T_l omega) = g(omega) a.e.), E(g(omega_1) - g(omega_0))^2 = 2E g^2 - 2E g Pg = 0 by invariance and harmonicity;
hence g(T_l omega) = g(omega) a.e. for every l with p > 0, in particular for l = L, which is the shift. So g is tau-invariant nu-a.e., hence beta-a.e.,
hence constant. The chain (tape, next label) is ergodic by the same argument as in 805. No positivity of 3-grams is needed, only (H2) and (P).
(3) *Counting.* As in 805 step 4: by ergodicity A_i/n, B_i/n -> 0 (zero drift (D) is exactly the statement E_nu[Delta a] = E_nu[Delta b] = 0, the
integrability being trivial since |Delta| <= 1) uniformly in i <= n, and I_n/n -> h a.s. Fix eps. The event E (|A_i| + |B_i| <= eps n for all i <= n and
I_n >= n(h-eps)) has probability >= 1/2 for large n, so there is one environment omega with P(E | omega_0 = omega) >= 1/2. Given omega, a label
sequence has probability exp(-I_n) and is applicable from omega when its probability is positive. So E contains >= (1/2) exp(n(h - eps)) sequences.
Let ell = ceil(eps n) + 2 + r and T_0 the 2 ell + 1 letters of omega around the head. By induction on i the finite tape T_0 after i moves is the part of
omega_i from ell + A_i letters left of the head to ell + B_i letters right of it, both numbers >= ell - eps n >= r + 2; every applicable label at omega_i
is applicable at the finite tape (L needs a left letter, R a right letter, M both neighbors), and the window of radius r that determines p lies inside it.
**The strip is not used here: T_0 is read as an unrestricted word.** Different label sequences from the same tape are different walks. Pad
T_0 to a tape T_1 of class C (first and last letter +, letter sum 1: put + c^i in front and c^i' + behind with one sign c) with <= 4 ell + 5 letters; the
sequences remain applicable (they read only the letters of T_0 and its neighbors inside T_0 on E). Lemma B (one build walk from the root to T_1,
<= 3 (4 ell + 5) moves) and Lemma R (return from the end tape, class C, <= 4 ell + 5 + eps n letters, <= 3 times that many moves) are statements about
the *unrestricted* machine (805 step 3; I reimplemented both and checked them on random class-C tapes up to 61 letters, FINAL section 3.6). The
resulting closed walks from the root are distinct (they differ in the middle part), have lengths N in [n, n + n'] with n' <= 30 eps n + 110, so for one
N in this range c_N >= exp(n(h - eps))/(2(n'+1)). With c_N <= mu^(N/2) this gives h <= (1/2) log mu after n -> infinity, eps -> 0. QED

**Remarks.** (i) 808 REPAIRS 3 asks the sweeps and the padding to stay inside the strip; this is unnecessary: only the *middle* part of the closed walk
is confined to the strip, and the counted object c_N is the unrestricted machine. This is the whole content of the "cheap closure" (L6).
(ii) The forbidden merge is part of the *machine* (an inapplicable label has no flow variable), so nu_J gives it probability zero by construction.
(iii) Factor 2: closed walks of length N correspond to order (N+2)/2, hence mu^(N/2): exp(2h) per move to growth per order. Confirmed on K = 3 (FINAL 3.1).
(iv) (D) is a *hypothesis on exact numbers*. An interval enclosing 0 (as in the certificates, +-6e-47) is a consistency check; exactness comes from the
algebraic construction (sections 3 and 4): exact rational parent plus exact dyadic multiples of exactly stationary cycle fields.

## 3. Proposition B (literal closed walks are stationary directions) [proved; same as 1211 Lemma 1 / 1008 / LONG_LITERAL, with the details]

A *configuration* x = (lev, word, head) is a finite word read from height lev with a marked head letter; it *fits* radius r if head <= r and
len - 1 - head <= r. A move x -> x' is *protected* (memory m) if lev, the first m letters, the last m letters and the letter sum are unchanged and x'
is in the strip (only the unprotected middle is rewritten). Let mu_x be the probability measure on Omega under which the letters around the head spell
x and the exterior follows beta conditionally: mu_x = beta(. | x) (a probability measure). By the Markov property the exterior law given x depends on
(lev, first m letters) on the left and (end level, last m letters) on the right, hence is the same for x and x' if the move is protected.
**Claim.** Let (x_e, l_e, x'_e)_e be finitely many protected edges with integer (or real) weights f(e), every *source* x_e fitting radius r, and
suppose the net flow vanishes at every configuration (sum over edges into x of f = sum over edges out of x). Put
J(w, l) := sum over edges e with l_e = l such that w completes x_e of f(e) beta(w)/beta(x_e).
Then J satisfies (S_r) and (D). For a closed literal walk with f = beta(anchor) the coefficient is beta(anchor)/beta(x_t) on the completions of x_t
(multiplicities add, a repeated edge counts twice).
**Proof.** (a) T_l maps the cylinder of x bijectively onto the cylinder of x' (protected, so the exterior letters are untouched), and by protection
(T_l)_# mu_x = mu_x'. (b) The signed measure nu_l := sum_{e: l_e = l} f(e) mu_{x_e} has density J(w,l)/beta(w) wrt beta: mu_x(C(w)) = beta(w)/beta(x) for a
window w completing x (x fits, so x is inside w), and given w the exterior is beta(. | w) because x is inside w. So nu_l = mu_J(., l). (c) sum_l (T_l)_# nu_l =
sum_e f(e) mu_{x'_e} and sum_l nu_l = sum_e f(e) mu_{x_e}; their difference is sum_x (in(x) - out(x)) mu_x = 0. By Proposition A this is (S_r) (r >= m).
(d) Drift: a(x) = head, b(x) = len - 1 - head are functions on configurations, each edge changes them by Delta(l_e), mu_x has mass 1 and the windows completing x_e carry
total J-mass f(e) (sum over completions of beta(w)/beta(x_e) = 1): sum J Delta = sum_e f(e)(a(x'_e) - a(x_e), b(x'_e) - b(x_e)) = sum_x (in(x) - out(x)) (a(x), b(x)) = 0. QED
Targets need not fit (they cancel against sources or against negative edges): this is the "signed diamond" class of 1211; **protection is essential**
(negative control: with only 2 or 1 protected letters at memory 3, 1 of 80 and 8 of 80 random 4-step cycles violate rows, worst residual 3e-3 and 0.46;
with the coefficient beta(w) instead of beta(w) beta(A)/beta(x_t), 20 of 80 violate, FINAL 3.4).

**Corollary B1 (context reduction, 1211 Prop. 2).** If every state of a cycle with base word A has the same first 3 and last 3 letters, then for a
letter ell, beta(ell x)/beta(x) = P(ell | first m letters) is the same for all states, so Gamma(cycle at A) = sum_ell Gamma(cycle at ell A) (the windows
completing x_t are partitioned by the extra letter). This is also the "children sum to the lifted parent" statement of the sparse note (L5); the child
cycles are themselves protected cycles (their first/last m letters are ell + the old first m-1), so each child is separately a stationary direction.

## 4. Proposition C (signed intermediate parent) [proved]

**Statement.** Let r0 >= m and let p be a *signed* flow at radius r0 satisfying (S_r0), (D) and sum p = 1. Let r >= r0, p_r its lift to radius r
(Prop. A), and c_1..c_k flows at radius r satisfying (S_r) and (D) (for example literal cycle fields, Prop. B). For real a_i put J = p_r + sum a_i c_i. Then J
satisfies (S_r) and (D) with sum J = 1 + sum a_i sum c_i. If J(w,l) > 0 for **every** legal (w,l) at radius r and M := sum J > 0 (automatic), then
J/M is a flow to which Theorem F' applies, and mu >= exp(2 h(J/M)) with h(J/M) = sum_{w,l} J log(nu_J/J) / M.
**Proof.** (S_r) and (D) are linear equations; the lift p_r satisfies them by Proposition A applied to the *signed* measure mu_p (the proof is linear in
J); cycle fields satisfy them by Prop. B. Positivity of the completed J is what Theorem F' needs; nothing is evaluated on p itself. h is computed from
the completed fine flows only. QED
**Where a coarse quantity could enter and does not.** (i) The entropy formula contains only J(w,l) and nu_J(w) = sum_l J(w,l) at the fine radius; no ratio
involves p. (ii) p_r(W,l) = p(center,l) beta(W)/beta(center) multiplies a signed rational by a *positive* cylinder ratio, so signs are those of p. (iii) The
only division by a possibly-negative number is in the producer, not in the argument. In an implementation, enclosing a signed rational numerator/denominator
needs endpoint order chosen by the sign of the numerator (the notes say this is done); that is an implementation point for the other referee (1253).

## 5. Proposition D (sparse refinement, inherited positivity, local entropy) [proved]

(a) *Lift preserves mass and entropy.* For F >= 0 stationary at radius r0 and F_r its lift, nu_r(W) = nu_0(c) beta(W)/beta(c) and F_r(W,l)/nu_r(W) = F(c,l)/nu_0(c)
(c = center), and sum over the completions W of c of beta(W)/beta(c) is 1, so S(F_r) = S(F) and H(F_r) = H(F), with S = sum F, H = sum_w nu log nu - sum_e F log F.
(b) *Local change.* Let delta be supported on a finite set Jset of fine flows, W the windows they touch, nu' = nu + delta-sums. Then
H(F+delta) - H(F) = sum_W [nu' log nu' - nu log nu] - sum_Jset [F' log F' - F log F], S' - S = sum_Jset delta. Both are sums over Jset and W only: obvious from the
definitions; the stable form h' - h = [sum_J delta (-log(F/nu) - h) - sum_J F phi(delta/F) + sum_W nu phi(delta_nu/nu)]/(S + delta S), phi(x) = (1+x) log(1+x) - x,
is an algebraic identity (checked to 3e-16, `t5b_sparse_identity.py`; end-to-end 4e-16 in `t9_sparse_chain.py`).
(c) *Inherited positivity.* If F is positive on all legal pairs at radius r0 then F_r > 0 on all legal pairs at r (the factor beta(W)/beta(c) > 0, and every legal
(W,l) has a legal center). For F_r + delta, positivity needs checking exactly on Jset, **provided Jset contains every fine flow where delta is nonzero**: this is
the one requirement the sparse verifier must satisfy (it must enumerate for every cycle edge all completions of the source state; a window missed by the
enumeration would silently inherit positivity it does not have). It is a property of the implementation, not of the mathematics.
(d) *Children.* Corollary B1 gives: the sum over all children of a parent cycle equals the lifted parent field exactly; any subset of children is a sum of stationary
directions. The exact sum identity was checked numerically to 4e-16 (`t5_signed.py`).
(e) *Chains of increments.* A finite sum of stationary fields is stationary; positivity of the final flow must be checked on the union of supports of all increments
(each window with the *cumulative* value, since the earlier increments change the value there). Entropy and mass are cumulative. The notes (1228 section 5) say the
checker reconstructs all prior increments inside every changed window; that is what this requires.
