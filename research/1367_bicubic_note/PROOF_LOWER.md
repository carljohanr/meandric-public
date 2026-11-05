# 1367: a lower bound for bicubic Hamiltonian cycles by stationary flows

**Result.** `mu^2 >= 10.086918840819` (radius-10 flow certificate `cert/flow_r10.bin`, checked by `check_flow_a.py`
and `check_flow_b.py`). Here `mu^2 = lim z_n^(1/n)`, `z_n` the number of bicubic planar maps with `2n` vertices, a
Hamiltonian cycle and a marked oriented cycle edge (602 Theorem A). With 602's upper bound:
`10.0869 <= mu^2 <= 10.1571`. Status: proved here (this file), not reviewed; the certificate is checked by two
programs by one author.

The method is Theorem F of 805 / Paper B Section 6, adapted to the two-stack machine. Three things differ from the
meander tape and make it simpler: the machine is *reversible* (its transfer matrix is symmetric), so detailed
balance gives exactly stationary flows with no linear solve; the ergodicity step uses the two push moves instead of a
head shift; and the closing step uses a charge invariant (Lemma S).

## 1. The machine (602, Sections 2 and 3)

State `x = (U, L, c)`: two finite stacks over {0,1} and a parity `c`. Moves: push A `(U c, L, 1-c)`, push B,
pop A `(U', L, 1-c)` when `U = U' y` with `y = 1-c`, pop B. `z_n = (A^(2n))[e, e]`, `e = (empty, empty, 0)`, and
`z_(m+n) >= z_m z_n`, so `mu^2 = sup_n z_n^(1/n)` and `z_n <= mu^(2n)` (Fekete).

**Lemma 1 (symmetry).** `y` is a successor of `x` iff `x` is a successor of `y`. Pushing `c` at parity `c` and
popping the top `c` at parity `1-c` undo each other. So `A` is symmetric, and a walk read backward is a walk.

**Lemma 2 (charge).** Let `N0`, `N1` count the zeros and ones in both stacks. Every move preserves
`N0 - N1 - c`: pushing `c` changes `N0 - N1` by `+1` if `c = 0` and `-1` if `c = 1`, which is the
change of `c`; popping `1-c` likewise. Hence every state reachable from `e` has `N0 - N1 = c`.

**Lemma S (sweep).** Every state `x = (U, L, c)` with `N0 - N1 = c` is joined to `e` by a walk of at most
`5|U| + 3|L|` moves. (With Lemma 1 also from `e` to `x`.)

*Proof.* A *transfer* from `U` to `L` (when the top of `U` is `c`) is push `c` on `L`, then pop `U`: two moves, the
top `c` moves from `U` to `L`, the parity is unchanged. Repeat, at parity `c`:
1. if `U` is nonempty with top `1-c`, pop `U`;  2. else if `L` is nonempty with top `1-c`, pop `L`;
3. else if `U` contains a symbol `1-c`, transfer from `U` to `L` (its top is `c`, by 1);
4. else transfer from `L` to `U`.
Rule 4 is possible unless `L` is empty: then `U` is `c^k`, and `N0 - N1 = c` forces `k = 0` and `c = 0`, so
`x = e`. The charge is preserved (Lemma 2), so this holds at every step. When rule 4 first applies, `U` consists of
the symbol `c` only; afterwards `U` only receives that symbol (by rule 4, at parities where all nonempty tops equal
the parity) and stays one-symbol, and at every later step where rules 1 and 2 fail the parity equals that symbol,
so rule 3 never applies again. So each original symbol of `U` is transferred at most once by rule 3, and each
symbol is transferred at most once by rule 4. Transfers do not change the total size, pops lower it by one, so
there are `|U| + |L|` pops and at most `|U| + (|U| + |L|)` transfers: at most `5|U| + 3|L|` moves. Each step
lowers `(size, |U|)` (rule 3) or `(size, |L|)` lexicographically, so the procedure ends, at `e`. ∎

`sweep_test.py` checks Lemmas 1, 2, S on all states of size at most 10 (5,193 states; worst ratio 2.6 moves per
symbol) and that breadth-first search reaches exactly the states of the charge rule.

## 2. Infinite stacks, windows and flows

`Omega`: pairs of one-sided infinite stacks `U = (U_1, U_2, ...)`, `L = (L_1, ...)` (index 1 = top) and a parity.
The same four moves act on `Omega` (pops need `U_1 = 1-c`). Fix a radius `r >= 1`. The *window* of `omega` is
`w = (U_1..U_r, L_1..L_r, c)`; the *extended window* keeps `r+1` symbols of each stack.

*Reference law.* `beta` is the law of a stationary Markov chain `(X_1, X_2, ...)` on {0,1}, read downward, of
memory `m <= r`, with all transition probabilities positive and *reversible* (its law is invariant under reversal of
finite blocks). The *continuation* of a word `w_1..w_k` (`k >= m`) is the conditional law of `X_(k+1), X_(k+2), ...`
given `X_1..X_k = w_1..w_k`. Here we use `m = 1`: `T(x, y) = s` for `x != y`, `1 - s` for `x = y`, `0 < s < 1`,
stationary law `(1/2, 1/2)`; a two-state chain is always reversible.

*Flow.* Numbers `J(w, xi) >= 0` for windows `w` and moves `xi` allowed at `w`, total `1`.
`nu(w) = sum_xi J(w, xi)`, `p(xi | w) = J(w, xi) / nu(w)`. The law `Pi` on `Omega`: the window has law `nu`, the two
far parts are independent continuations of the window's stacks. The process: `omega_0 ~ Pi`, then move `xi_i` with
probability `p(xi_i | w(omega_i))`, `omega_(i+1) = xi_i(omega_i)`.

Conditions:
- **(S)** for every extended window `Y`: `sum_(w, xi) J(w, xi) P_(w,xi)(Y) = nu(Y°) beta(Y_(U,r+1) | Y°_U) beta(Y_(L,r+1) | Y°_L)`,
  where `P_(w, xi)(Y)` is the probability, with the far parts continued from `w`, that the move gives extended
  window `Y`, and `Y°` is `Y` without the symbols at depth `r+1`.
- **(D)** total flow of push A = total flow of pop A; the same for B.
- **(P)** `J(w, push A) > 0` and `J(w, push B) > 0` for every window `w`.

**Theorem L.** Under (S), (D), (P): `mu >= e^h` with `h = sum_(w,xi) J log(nu(w)/J(w,xi))`, so `mu^2 >= e^(2h)`.

*Proof.* **Step 1 (stationarity).** Let `omega ~ Pi`, `xi ~ p(. | w)`, `omega' = xi(omega)`. Conditionally on `w`,
`xi` and the extended window `Y` of `omega'`, the symbols of `omega'` below depth `r+1` are the continuations of
`Y`'s stacks: after a push on `U` they are the old symbols `U_(r+1), ...`, which continue `U_1..U_r`, the last `m`
symbols of `Y`'s `U`-part (here `m <= r` is used); after a pop on `U` the symbols `U_(r+1), U_(r+2)` became part
of `Y` and the rest continues them; the untouched stack keeps its continuation beyond its revealed symbol at depth
`r+1`. By (S) the law of `Y` is `nu(Y°)` times the continuation probabilities of the two depth-`(r+1)` symbols,
which is the law of the extended window under `Pi`. So `omega' ~ Pi`.

**Step 2 (ergodicity).** By (P), `nu > 0`, so `Pi` and `P := beta x beta x (1/2, 1/2)` have the same null sets; the
moves map `P`-null sets to `P`-null sets. Let `f` be bounded with `E[f(omega_1) | omega_0] = f(omega_0)`. Stationarity
gives `E(f(omega_1) - f(omega_0))^2 = 0`, so `f(xi omega) = f(omega)` a.s. for every move of positive probability, in
particular for push A and push B at every window. Write `g(U, L) = f(U, L, 0)`. Push A and push B at parity 0 and 1
give, a.s., `f(U,L,0) = f(0U,L,1) = f(U,0L,1)` and `f(U,L,1) = f(1U,L,0) = f(U,1L,0)`. Hence
`g(1U, L) = g(U, 1L)` and `g(0U, L) = f(0U, 0L, 1) = g(U, 0L)`: a.s. `g(xU, L) = g(U, xL)` for `x` in {0,1}.
Glue the stacks into one two-sided sequence `Z = (..., U_2, U_1 | L_1, L_2, ...)`; the identity says that `g`, as a
function of `Z`, is invariant under the shift of `Z` by one place. Under `P` the two halves of `Z` are independent
copies of the chain read away from the cut. Because `beta` is reversible with positive transitions, this law and the
stationary two-sided chain have the same null sets (they differ only in the dependence of `m` symbols across the
cut). The stationary irreducible chain is ergodic under the shift, so `g` is a.s. constant, and so is
`f(U, L, 1) = g(1U, L)`. A stationary Markov process whose bounded harmonic functions are constant is ergodic.

**Step 3 (counting).** Let `H^A_i, H^B_i` be the height changes of the two stacks after `i` moves and
`I_t = sum_(i<t) -log p(xi_i | w(omega_i))`. By ergodicity and (D), `H_t/t -> 0`, so `max_(i<=t) |H^A_i| + |H^B_i| <= eps t`
eventually, and `I_t/t -> h`. For large `t` this event `G_t` (with `I_t >= t(h - eps)`) has probability at least
`1/2`, hence conditional probability at least `1/2` given `omega_0` for some fixed `omega_0`. Given `omega_0` the moves
determine the path and a move sequence has probability `exp(-I_t)`, so at least `(1/2) e^(t(h-eps))` move sequences
lie in `G_t`. Let `D = ceil(eps t) + 1`, `alpha`, `gamma` the top `D` symbols of the two stacks of `omega_0`, and
`x = (P alpha, gamma, c_0)` with `P` a run of zeros or ones of length at most `2D + 1` chosen so that
`N0 - N1 = c_0` (Lemma 2). Every sequence in `G_t` is a legal walk of the finite machine from `x`, since it never
looks below depth `eps t + 1`. By Lemma S and Lemma 1 there is a walk `e -> x` of length at most `5|x| = O(eps t)`,
and every end state `y` (size at most `|x| + 2 eps t`, and reachable) closes to `e` in `O(eps t)` moves. The
resulting closed walks are distinct (the prefix is fixed) and have lengths in `[t, t + C eps t + C']`, so for some
`n` in that range of half-lengths `z_n >= e^(t(h-eps)) / (2(C eps t + C' + 1))`. With `z_n <= mu^(2n)`,
`log mu >= (h - eps) t / (t + C eps t + C') - o(1)`; let `t -> infinity`, then `eps -> 0`. ∎

## 3. Reversible flows

**Proposition R.** Let `beta` be as above (memory `m <= r`, here `m = 1`), and let `X_A(w) > 0` be given for the
windows with `a_1 != c` and `X_B(w) > 0` for those with `b_1 != c`. Put `beta(w) = beta(a) beta(b) / 2`, with
`beta(a)` the stationary probability of the word `a` and `q(x | a) = P(X_0 = x | X_1..X_m = a_1..a_m)` (for the
symmetric chain, `q(x | a) = T(a_1, x)`), and
- `J(w, pop A) = beta(w) X_A(w)`,  `J(w, push A) = beta(w) q(c | a) X_A(w^A)`, where `w^A = ((c, a_1..a_(r-1)), b, 1-c)`,
- the same for B. Then `J` (normalized) satisfies (S), (D) and (P).

*Proof.* On `Omega`, `Pi` has density `rho(omega) = nu(w)/beta(w)` with respect to `P`. Push A maps
`omega = (U, L, c)` to `omega' = (cU, L, 1-c)`, and `P(d omega') = q(c | U_1..U_m) P(d omega)` under this map
(stationary chain: adding a symbol on top multiplies the probability by the upward transition, which equals the
downward one by reversibility). Pop A at `omega'` is the inverse map. The flow through the edge is
`rho(omega) p(push A | w) P(d omega) = beta(w)^(-1) J(w, push A) P(d omega) = q(c | a) X_A(w^A) P(d omega)` one way
and `beta(w^A)^(-1) J(w^A, pop A) P(d omega') = X_A(w^A) q(c | a) P(d omega)` the other way: equal. So `Pi` with
these moves is in detailed balance; inflow equals outflow at every configuration, which is (S) (integrate over the
configurations with extended window `Y`). Summing the edge identity gives (D). (P) holds because all `X > 0`. ∎

So any positive `X` gives a bound; optimization only chooses `X`. The exact flows are integers after scaling:
with `s = sn/sd`, multiply by `8 sd^(2r-1)` (checker A) or `8 sd^(2r+1)` (checker B).

## 4. The certificates

Reference law: symmetric chain, `s = 41/50` (memory 1; the stacks prefer to alternate). `X` found by the fixed-point
iteration `fixpoint.py` (at an optimum, `log X(v)` is the weighted mean of `log(nu/coefficient)` over the edges that
use `X(v)`), scaled so the largest is `2^60`, rounded, at least 1. File format in `make_cert.py`.

| r | file | sha256 (first 16) | windows | `e^(2h) >=` (both checkers) |
|---|---|---|---|---|
| 4 | `cert/flow_r4.bin` | c4023c79ad842835 | 512 | 10.030980878326 |
| 8 | `cert/flow_r8.bin` | b699595c121609ec | 131,072 | 10.079276064125 |
| 10 | `cert/flow_r10.bin` | 20bcd572c19bb5c6 | 2,097,152 | **10.086918840819** |

- Checker A (`check_flow_a.py`): builds `J` from Proposition R in exact integers, checks (P) and (D) exactly,
  encloses `h` in 160-bit interval arithmetic (mpmath `iv`).
- Checker B (`check_flow_b.py`, separately written, stacks as strings, own reader): checks (S) on all
  `2^(2r+3)` extended windows (8,388,608 at `r = 10`) in exact integers without using Proposition R, checks (D)
  and (P), and bounds `h` with Python `decimal` (correctly rounded `ln`, one unit in the last place added or
  subtracted, sums rounded toward the safe side).
- Teeth (`logs/teeth_flow_b.log`, radius 8): omitting the reference factor in the pushes, changing one push value
  by one unit, or revealing one symbol instead of two after a pop: (S) fails (524,288, 6 and 524,288 unbalanced
  windows).
- Float values by radius (`logs/float_table.log`, `s = 41/50`): 9.497, 9.880, 9.988, 10.031, 10.053, 10.065,
  10.073, 10.079, 10.084, 10.087 for `r = 1..10`; with the iid reference law `s = 1/2`: 9.135, 9.577, 9.735, 9.824,
  9.878, 9.917, 9.944, 9.966, 9.983 for `r = 1..9`. The gap to the estimate 10.113 shrinks slowly (0.026 at `r = 10`).
