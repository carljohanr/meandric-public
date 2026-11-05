# 602: why a verified row certificate bounds the number of Hamiltonian cycles on bicubic planar maps

This note proves three statements and describes the certificate that the verifiers check.

- **Theorem A (bijection).** Rooted bicubic planar maps with `2n` vertices and a Hamiltonian cycle are in
  bijection with the words of length `2n` of the two-stack language `W_n` (Section 2).
- **Theorem B (transfer).** `|W_n| = (A^(2n))[e, e]` for the transfer `A` of Section 3, `e` the empty state.
- **Theorem C (bound).** If a certificate `(r, a, p, t, P/Q, u)` passes the verifier, then
  `z_n <= (P/Q)^(2n)` for every `n >= 1`, hence `mu^2 := lim z_n^(1/n) <= (P/Q)^2`.

Status labels: **[proved here]**, **[relied on]**, **[test]** (evidence, not part of the proof).
It is written against the verifiers `verify_rows.py` and `verify_slow.cpp` (string views), not against the solver.

## 1. The objects

A *planar map* is a connected finite multigraph (loops and multiple edges allowed) embedded in the oriented
sphere, considered up to orientation-preserving homeomorphism. It is *cubic* if every vertex has degree 3 and
*bicubic* if it is cubic and bipartite. A *Hamiltonian cycle* is a cycle of the map through every vertex exactly
once. A *rooted configuration* of size `n` is a triple (bicubic planar map with `2n` vertices, Hamiltonian
cycle `H`, one edge of `H` marked and oriented). `z_n` is the number of rooted configurations of size `n`.

Remarks.
- A bicubic map has as many white as black vertices (count the edges from each color class), so the number of
  vertices is even.
- A connected bipartite map has exactly two proper colorings, exchanged by swapping the colors. So a rooted
  configuration is the same thing as a *bicolored* map with a Hamiltonian cycle and a marked (unoriented) edge
  of the cycle: orient the marked edge from its black end to its white end. This is the convention of the
  literature ("planar bicubic maps with 2N vertices endowed with a Hamiltonian cycle and a marked visited edge",
  Di Francesco, Duplantier, Golinelli, Guitter, arXiv:2210.08887). **[relied on]** for the identification of
  `z_n` with the published numbers; **[test]** the first ten values computed here from the transfer,
  `2, 8, 40, 228, 1424, 9520, 67064, 492292, 3735112, 29114128`, are the published ones (OEIS A116456, Table 1
  of Guitter, Kristjansen, Nielsen 1999).
- Rooted configurations have no nontrivial automorphism. The number of unrooted configurations lies between
  `z_n / (4n)` and `z_n`, and the automorphism-weighted sum of the literature is `z_n / (4n)` (GKN eq. 4.5:
  `M_N^E = 2 * 2N * Z_N`). None of this changes the growth rate.

## 2. Words  [proved here]

Letters: `A(`, `A)`, `B(`, `B)`. For a word `w = w_0 ... w_(N-1)` the *A-subword* is the word over `(`, `)`
read off its A-letters, likewise the B-subword. `W_n` is the set of words of length `2n` such that

1. the A-subword and the B-subword are Dyck words, and
2. for every matched pair of parentheses (in the A-subword or in the B-subword), the two positions `i < j` of
   its letters in `w` have opposite parity (`j - i` is odd).

**Theorem A.** `z_n = |W_n|`.

*Proof.* Let `(m, H, e)` be a rooted configuration. Number the vertices `0, 1, ..., 2n-1` along `H`, starting
at the head of the root edge `e` and following its direction; so `e` runs from vertex `2n-1` to vertex `0`.

(i) *The third edges.* `H` is a simple closed curve on the sphere; its complement consists of two open disks
(Jordan, Schoenflies), the one on the left of `H`, called `D_A`, and the one on the right, `D_B`. Every vertex
has degree 3 and two of its half-edges are on `H`, so exactly one half-edge is not. An edge outside `H` cannot
be a loop (a loop uses two half-edges of one vertex outside `H`). Hence the `n` edges outside `H` form a perfect
matching of the vertices, and each of them is an arc inside `D_A` or inside `D_B` (*chord*).

(ii) *The word.* Give vertex `j` the letter `A(` or `A)` if its chord lies in `D_A` and `j` is its smaller or
larger end, and `B(`, `B)` likewise for `D_B`. Two chords in the same disk are disjoint arcs joining boundary
points of a disk, so their ends do not interlace (`i < k < j < l` is impossible for chords `{i, j}`, `{k, l}`).
A non-crossing perfect matching of a linearly ordered set, with smaller ends written `(` and larger ends `)`,
is a Dyck word together with its parenthesis matching. This is condition 1.

(iii) *Bipartite = condition 2.* Color vertex `j` by the parity of `j`. The edges of `H` join `j` to `j+1` and
`2n-1` to `0`, always opposite parities because `2n` is even. A proper coloring of the connected map restricts
to a proper coloring of the cycle `H`, which is the parity coloring or its swap. So the map is bipartite iff
the parity coloring is proper iff every chord joins two positions of opposite parity.

So the configuration gives a word `Phi(m, H, e)` in `W_n`.

(iv) *Inverse.* Given `w` in `W_n`, put `2n` points on a circle oriented counterclockwise, draw the matched
A-pairs as non-crossing chords inside and the matched B-pairs as non-crossing arcs outside, and root the circle
at the edge from point `2n-1` to point `0`. This is a connected cubic planar map with a Hamiltonian cycle; by
(iii) it is bipartite; its word is `w`.

(v) *Injectivity.* A planar map is determined, up to orientation-preserving homeomorphism, by its rotation
system (the cyclic order of the half-edges around each vertex). In a rooted configuration the rotation at vertex
`j` is (half-edge of `H` toward `j-1`, half-edge of `H` toward `j+1`, third half-edge) in clockwise or
counterclockwise order according to whether the third half-edge lies in `D_A` or `D_B`, and the third
half-edges are paired by the parenthesis matchings. All of this is read off the word, so two configurations
with the same word are equal. ∎

**Lemma 2.1 (growth rate).** `z_(m+n) >= z_m z_n`, so `mu^2 := lim z_n^(1/n) = sup_n z_n^(1/n)` exists; and
`z_n <= 16^n`.

*Proof.* The concatenation of `w` in `W_m` and `w'` in `W_n` is in `W_(m+n)`: the positions of `w'` are shifted
by the even number `2m`, so parities are preserved, and matchings do not mix. Concatenation is injective because
the lengths are fixed. Fekete's lemma gives the limit. `W_n` is a subset of the `4^(2n)` words. ∎

## 3. The transfer  [proved here]

A *stack* is a word over the parities `0`, `1` (bottom first, top last). A *state* is a triple `(U, L, c)` of
two stacks and a parity `c`. The moves from `(U, L, c)` are

- **oA**: `(U c, L, 1-c)` (push `c` on `U`); always available;
- **oB**: `(U, L c, 1-c)`; always available;
- **cA**: `(U', L, 1-c)` where `U = U' x`; available iff `U` is not empty and its top `x` is `1-c`;
- **cB**: `(U, L', 1-c)` where `L = L' x`; available iff `L` is not empty and its top `x` is `1-c`.

The four moves change `(|U|, |L|)` by `(+1, 0)`, `(0, +1)`, `(-1, 0)`, `(0, -1)`, so the successors of a state
are distinct. `A[x, y] = 1` if `y` is a successor of `x`, else `0`; `(A^N)[x, y]` is the number of walks of `N`
moves from `x` to `y`. `e = (empty, empty, 0)`.

**Theorem B.** `(A^(2n))[e, e] = |W_n|`, and `(A^N)[e, e] = 0` for odd `N`.

*Proof.* Read the moves oA, oB, cA, cB as the letters `A(`, `B(`, `A)`, `B)`. A walk from `e` is determined by
its sequence of moves, and after `j` moves the parity component is `j mod 2`, the parity of the position of
the next letter. By induction on `j`: after the first `j` letters, `U` lists, from bottom to top, the parities
of the positions of the unmatched `A(` in increasing order of position (and `L` the same for `B`). Indeed
`A(` at position `j` appends the parity `j mod 2 = c`. A letter `A)` at position `j` is legal in a Dyck prefix
iff there is an unmatched `A(`; its partner in the parenthesis matching is the last unmatched `A(`, which is the
top of `U`; condition 2 for this pair says that the top differs from `j mod 2 = c`. So the move cA is available
exactly when the letter `A)` keeps the A-subword a Dyck prefix and satisfies condition 2, and then it removes
the top. Hence the sequences of `N` moves that are walks from `e` are exactly the words of length `N` whose A-
and B-subwords are Dyck prefixes with all closed pairs satisfying condition 2; the walk ends with both stacks
empty iff both subwords are Dyck words; then `N` is even and the final parity is `0`. ∎

So `z_n = (A^(2n))[e, e]`.

*Remark (not used).* On states reachable from `e` the number of `0`s minus the number of `1`s in the two stacks
together equals `c` (every move at parity `c` changes this number by `(-1)^c`). The bound below is proved for
all states, reachable or not; discarding the unreachable rows did not change the float value at radius 3 and 4.

## 4. Sub-invariant positive functions  [proved here]

**Lemma 4.1.** Let `v > 0` be a function on all states and `R > 0` with
`(A v)(x) := sum over the successors y of x of v(y) <= R v(x)` for every state `x`. Then
`(A^N)[x, x] <= R^N` for all `x` and `N >= 0`.

*Proof.* Every state has at most four successors, so all sums are finite. By induction
`(A^N v)(x) <= R^N v(x)`, using `A >= 0`. Then `(A^N)[x, x] v(x) <= (A^N v)(x) <= R^N v(x)`; divide by
`v(x) > 0`. ∎

With `x = e`: `z_n <= R^(2n)` for every `n`, so `mu^2 <= R^2`. The bound is per size; no asymptotic input.

## 5. The weight and the rows  [proved here]

Fix a radius `r >= 2`, positive rationals `a, p, t`, and a table `u` of positive integers indexed by contexts.

- For a stack `S`: `E(S)` = number of `i` with `S_i = S_(i+1)` (adjacent equal pairs), `T(S)` = number of `i`
  with `S_i = S_(i+1) = S_(i+2)` (adjacent equal triples).
- `view_k(S)` = the string of the last `k` symbols of `"^" + S` (the whole string if it is shorter). So a view
  is either `k` parities, or the bottom mark `^` followed by fewer than `k` parities (then it is the whole
  stack). `V_k` denotes the set of all such strings; it has `2^(k+1) - 1` elements, and every element of `V_k`
  is the view of some stack.
- The *context* of `(U, L, c)` is `(view_r(U), view_r(L), c)`; the *extended views* are `EU = view_(r+1)(U)`,
  `EL = view_(r+1)(L)`. Equivalently `EU` is the view with one more symbol revealed under it (`^`, `0` or `1`)
  when the view shows `r` parities, and `EU = view_r(U)` when the view shows the bottom mark.
- `v(U, L, c) = a^(|U|+|L|) p^(E(U)+E(L)) t^(T(U)+T(L)) u(view_r(U), view_r(L), c)`.

The table is indexed by all of `V_r x V_r x {0, 1}`, so `v` is defined and positive on all states. A *row* is
an element `(EU, EL, c)` of `V_(r+1) x V_(r+1) x {0, 1}`; every row is the triple (extended views, parity) of
some state, and conversely. There are `2 (2^(r+2) - 1)^2` rows.

For a string `s` in `V_(r+1)` (or such a string with one parity appended) let `view(s)` be `s` if it has at
most `r` symbols and the string of its last `r` symbols otherwise; let `pairs(s)` and `triples(s)` be the
numbers of adjacent equal pairs and triples of parities in `s` (the mark `^` is never part of one).

**Row lemma.** Let `(U, L, c)` be a state with extended views `EU`, `EL`.

1. *oA.* The new stack is `U c`. Its view is `view(EU c)`. `E(U c) - E(U) = pairs(EU c) - pairs(EU)` and
   `T(U c) - T(U) = triples(EU c) - triples(EU)`; both are 0 or 1.
2. *cA.* Available iff the last symbol of `EU` is a parity different from `c`. Let `EU-` be `EU` without its
   last symbol. The new stack `U'` has view `EU-`, and `E(U') - E(U) = pairs(EU-) - pairs(EU)`,
   `T(U') - T(U) = triples(EU-) - triples(EU)`; both are 0 or -1.
3. *oB, cB.* The same with `L`, `EL`.
4. In each move the other stack, its view and its numbers `E`, `T` are unchanged, `c` becomes `1-c`, and
   `|U| + |L|` changes by `+1` (oA, oB) or `-1` (cA, cB).

*Proof.* `EU` is the string of the last `r+1` symbols of `"^" + U`, or all of it.
(1) The last `r` symbols of `"^" + U c` are the last `r` symbols of `EU c`, because `EU` contains the last
`r - 1` symbols of `"^" + U` (or all of it); if `"^" + U c` has at most `r` symbols then `EU c` is that whole
string. The pairs of `U c` are those of `U` and the pair formed by the top of `U` and `c` (if `U` is not empty);
the triples of `U c` are those of `U` and the one formed by the two top symbols of `U` and `c` (if `|U| >= 2`).
The two last symbols of `"^" + U` are in `EU` because `r + 1 >= 2`; a `^` among them means that the pair or
triple does not exist, and `pairs`, `triples` do not count it.
(2) The top of `U` is the last symbol of `"^" + U`, hence the last symbol of `EU`; it is `^` iff `U` is empty.
`"^" + U'` is `"^" + U` without its last symbol, so its last `r` symbols (or all of it) form `EU-`, which
has at most `r` symbols. The pairs of `U` are those of `U'` and the pair of the two top symbols of `U` (if
`|U| >= 2`); the triples of `U` are those of `U'` and the one of the three top symbols (if `|U| >= 3`). The
three last symbols of `"^" + U` are in `EU` because `r + 1 >= 3`.
(3), (4) are clear. ∎

Consequently, for a state `x = (U, L, c)` with row `(EU, EL, c)`,

    (A v)(x) / ( a^(|U|+|L|) p^(E(U)+E(L)) t^(T(U)+T(L)) )
        =   a p^dE1 t^dT1 u(view(EU c), view(EL), 1-c)  +  a p^dE2 t^dT2 u(view(EU), view(EL c), 1-c)
          + [cA] (1/a) p^dE3 t^dT3 u(EU-, view(EL), 1-c)  +  [cB] (1/a) p^dE4 t^dT4 u(view(EU), EL-, 1-c),      (*)

where the exponents are the differences of `pairs` and `triples` named in the Row lemma. The right side is a
function of the row alone; call it `LHS(EU, EL, c)`.

**Covering lemma.** `A v <= R v` holds at every state if and only if `LHS(EU, EL, c) <= R u(view(EU),
view(EL), c)` holds for every row.

*Proof.* Divide `(A v)(x) <= R v(x)` by the positive number `a^(...) p^(...) t^(...)` and use (*). Every state
has a row and every row comes from a state. ∎

Compared with the meander transfer (262) there is nothing to bound over unseen symbols: a move changes only
the top of one stack, the revealed symbol is part of the row, and the counted patterns have length at most 3.
The three choices of the revealed symbol under each view are separate rows, so in effect the pop term enters
with the largest of its (up to) three destinations.

**Theorem C.** If every row satisfies `LHS <= (P/Q) u(context)`, then `z_n <= (P/Q)^(2n)` for all `n >= 1` and
`mu^2 <= (P/Q)^2`.

*Proof.* Covering lemma, Lemma 4.1 with `x = e`, Theorem B, Theorem A. ∎

## 6. Exact arithmetic and the certificate file

`a = an/ad`, `p = pn/pd`, `t = tn/td`. Multiply the row inequality by `an ad pn pd tn td Q`. A term
`a^s p^e t^f u` (`s = +-1`; `e, f` in `{-1, 0, 1}`) becomes `A(s) B(e) C(f) u Q` with `A(1) = an^2`,
`A(-1) = ad^2`, `B(1) = pn^2`, `B(0) = pn pd`, `B(-1) = pd^2`, and `C` likewise from `tn, td`; the right side
becomes `P an ad pn pd tn td u(context)`. `verify_slow.cpp` and `bic check` evaluate both sides in unsigned
128-bit integers after a bound that excludes overflow; `verify_rows.py` uses Python `Fraction`s on the unscaled
form.

File format `BICCERT1`: 8 bytes magic, `int32 r`, then `uint64 an, ad, pn, pd, tn, td, P, Q, count` with
`count = 4^(r+1)`, then `count` little-endian `uint32`. A view is coded by the integer whose binary digits are
a leading 1 followed by the displayed parities (top = lowest bit): `"^"` is 1, `"^01"` is `0b101`, a view of `r`
parities `w` is `2^r + w`. The entry at position `code(X) 2^(r+1) + code(Y)` is `u(X, Y, 0)`; positions with a
code 0 are unused. The entries for parity 1 are not stored: the reader defines

    u(X, Y, 1) := u(swap X, swap Y, 0),      swap = exchange of the parities 0 and 1 in the string.

This is only an encoding of a full table `u` on `V_r x V_r x {0, 1}`: the verifiers check every row of the
expanded table, and nothing in Sections 1-5 uses a symmetry. (The transfer commutes with the exchange of the
parities together with `c -> 1-c`, and with the exchange of the two stacks; the solver searches among tables
with both symmetries, which costs nothing in the float values.) The verifiers require that the number of
entries is `4^(r+1)`, that all used entries are positive and that the file has no trailing bytes.

## 7. Two remarks on the choice of the weight

**The exponents of `p` and `t` are not seen by the table.** `E(U)` is not a function of the context, but its
change under a move is a function of the row. That is why `p` and `t` enlarge the family, like `a` itself.

**Parameters that are odd under a symmetry cannot help.  [proved here, short]** Write the parameters as
`theta` (logarithms) and the table as `g = log u`. Each row inequality reads
`log sum_i exp(l_i(theta, g)) - g(context) <= log R` with `l_i` affine, so
`F(theta, g) = max over rows of the left side` is jointly convex and `rho(theta) = inf_g F(theta, g)` is convex.
For the factor `d^(|U| - |L|)` the exchange of the two stacks gives `rho(log d) = rho(-log d)`, hence
`rho(0) <= rho(log d)`: `d = 1` is optimal. The same argument applies to separate strand factors for the two
parities (exchange of the parities).

## 8. The certificates and the independent evidence

**Certified values.** (`cert/`, logs in `logs/cert_r*.log`, hashes in `cert/SHA256SUMS`.)

| radius | a | p | t | P/Q | (P/Q)^2 exactly | rows | checked by |
|---|---|---|---|---|---|---|---|
| 12 | 1783/2000 | 137/200 | 97/125 | 3187011595/10^9 | 10.157042906664444025 | 536,805,378 | `bic check`, `verify_slow` |
| 7 | 1771/2000 | 707/1000 | 773/1000 | 3191740268/10^9 | 10.187205938372711824 | 522,242 | `bic check`, `verify_slow`, `verify_rows.py` |

So, by Theorem C, `z_n <= (3187011595/10^9)^(2n)` for every `n >= 1` and `mu^2 <= 10.157042906664444025 < 10.15705`.
`P` is the smallest integer for which all rows hold with `Q = 10^9` (lowering it by one is rejected).

**Tests.**  [test]

- *Three row generators agree.* `bic check` (integer view codes, formulas for the pattern changes),
  `verify_rows.py` and `verify_slow.cpp` (strings, literal counting) produce the same rows: the row digest, a
  sum over rows of a hash of (context index, the destination indices, the exponent changes), is
  `7fa9757044f0d59d` for all three at radius 7 and `9b544869d96c399d` for the two C++ programs at radius 12.
- *Exact test on complete states.* `falsify_exact` applies the literal transfer to complete states, recounts
  `E` and `T` on the complete stacks before and after each move, reads the context off the complete stacks and
  compares in exact integers. Radius 12: all 83,886,082 states with `|U| + |L| <= 20` and 30,000,000 random
  states (long stacks of five kinds, one short stack, both heights within 3 of the radius): no violation. The
  largest ratio found on states, 3.187011594348 at `U = 0`, `L = 1111000000`, `c = 0`, equals the largest row
  ratio, so the bound of the table is attained by a real state. The program also recomputes `z_1..z_10` with
  the same successor routine.
- *Teeth.* Damaged certificates (P lowered by one unit; one entry lowered by 0.1 percent; two different entries
  raised by 0.1 percent) are rejected by every row checker (`logs/cert_r*.log`). The test on complete states
  (run here only up to 14 arches) rejects all four at radius 7 and three of the four at radius 12; the fourth
  damage sits in a context with two views of 12 parities, which no state with fewer than 24 arches has.
- *Per-size check.* `z_n <= (P/Q)^(2n)` holds for the 34 published values (`check_counts.py`).
