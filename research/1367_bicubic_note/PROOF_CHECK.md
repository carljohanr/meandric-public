# 1367: why `check_states` verifies the radius-12 certificate of 602

Objects as in 602 `PROOF.md` (refereed in 602r): states `x = (U, L, c)` (two stacks over {0,1}, written bottom
first, and a parity), the transfer `A` with moves push A (`U -> U c`), push B, pop A (allowed iff `U` is nonempty
and its top differs from `c`), pop B, each flipping `c`; the weight
`v(x) = a^(|U|+|L|) p^(E(U)+E(L)) t^(T(U)+T(L)) u(view_r(U), view_r(L), c)`, where `E` counts adjacent equal
pairs, `T` adjacent equal triples, and `view_r(S)` is the string of the last `r` symbols of `"^" + S`.
602 Lemma 4.1 and Theorems A, B: if `(A v)(x) <= R v(x)` at every state, then `z_n <= R^(2n)` for all `n`.

**Lemma T (truncation).** Let `r >= 2`. For a state `x = (U, L, c)` let `x'` keep only the top `r+1` symbols of
each stack (all of a stack that is shorter). Then `(A v)(x) / v(x) = (A v)(x') / v(x')`.

*Proof.* The moves allowed at `x` and `x'` are the same, because legality depends only on emptiness and the top
symbol, and a stack is empty iff its truncation is. Match the successors `y` of `x` with the successors `y'` of
`x'` by move. For each move, `v(y)/v(x) = a^(ds) p^(dE) t^(dT) u(y)/u(x)`, so it suffices that `ds`, `dE`, `dT`,
`u(x)` and `u(y)` agree.
- `ds = +-1` for both.
- A push of `c` on `U` adds a pair iff the old top equals `c`, and a triple iff the two old top symbols equal `c`
  (and exist). A pop removes the pair of the two top symbols iff they are equal, and the triple of the three top
  symbols iff they are equal. These involve at most the three top symbols, which `U` and its truncation share
  because `r + 1 >= 3`; existence (length at least 2 or 3) also agrees, since the truncation has length
  `min(|U|, r+1)`.
- `view_r(U) = view_r(U')`: if `|U| <= r+1` then `U' = U`; otherwise both show the top `r` parities.
- After a push the new view is the pushed symbol and the top `r-1` symbols; after a pop it is symbols `2..r+1`
  from the top (or the whole remaining stack with the mark). Both lie in the top `r+1` symbols. The other stack
  and its view are unchanged. ∎

So `A v <= R v` holds at all states iff it holds at the `2 (2^(r+2) - 1)^2` states with `|U|, |L| <= r+1`
(at `r = 12`: 536,805,378 states, the same number as the rows of 602). `check_states.c` enumerates exactly these
states, computes successors by the literal moves, recounts `E` and `T` on the whole stacks of `x` and of each `y`,
reads all views off the stacks, and compares `sum_y a^ds p^dE t^dT u(y) <= (P/Q) u(x)` after multiplication by
`an ad pn pd tn td Q`, in unsigned 128-bit integers; a start-up bound (largest entry times the largest coefficients,
times 4) shows that no product or sum exceeds `2^127`. Entries for parity 1 are read as
`u(X, Y, 1) = u(swap X, swap Y, 0)` (602 file format); every used entry must be positive.

It shares no code with `bic.cpp`, `verify_slow.cpp`, `verify_rows.py` or `falsify_exact.cpp` of 602, which were not
read. It shares the file format and the mathematical definitions, which come from 602 `PROOF.md` Sections 3 to 6.
