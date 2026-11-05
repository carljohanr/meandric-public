# Linear-record completion of a bounded-winding half

**Theorem.** Let a connected full meander be cut at a gap, and let a resulting half-state have `w` unmatched through-arch stubs in total. If its ordinary vertices number `2m`, then the number `H_m(i,j)` of such half-states with `i` upper and `j` lower stubs (`i+j=w`) obeys

`H_m(i,j) <= 4^(w-1) C^(3w/2) M_m`,

where `C=sup_n M_(n+1)/M_n<infinity`. Consequently, uniformly in `m,n,w>=2`,

`#{order-(m+n) meanders with at most w arches crossing the cut} <= exp(Bw) M_m M_n`

for an absolute `B`. This improves the `exp(O(w log(w+1)))` inverse-record loss of the arbitrary-junction completion in `1340_alpha_campaign/PROOF_BOUNDED_WINDING_COMPLETION.md`.

**Initial closure.** Use the same canonical closure as that proof: append `w` suffix vertices, connect upper stubs in reverse order to the first `i`, lower stubs in reverse order to the final `j`, and fill opposite-color suffix edges by adjacent short arches. The resulting system is planar, and every loop contains at least one of the `w` original suffix vertices. Label these vertices `v_1,...,v_w` in river order. They are all consecutive before mergers.

**Fixed left-to-right merger schedule.** Maintain the processed prefix of the suffix through original vertex `v_k`, together with all wiggle vertices inserted at earlier stages. Maintain the stronger invariant that all processed-prefix vertices lie on one loop and that its rightmost vertex `r_k` is adjacent to the next untouched original suffix vertex `v_(k+1)`. Initially `r_1=v_1`. At stage `k=1,...,w-1`, inspect the loops of `r_k` and `v_(k+1)`.

* If they are already on the same loop, record `skip` and set `r_(k+1)=v_(k+1)`.
* Otherwise use the unique one of the CO, OO, CC junction cases from `113_ratio_bound/PROOF.md`, Lemma 5, at this adjacent pair. Record that type. The move merges exactly those loops. In CO the next original vertex stays untouched and is the new rightmost representative. In OO, the mirror wiggle replaces `r_k` by three vertices immediately before the untouched `v_(k+1)`; the subsequent CO exchange uses the third replacement and `v_(k+1)`, so again `r_(k+1)=v_(k+1)`. In CC, the ordinary wiggle replaces `v_(k+1)` by three vertices; the CO exchange uses the old `r_k` and the first replacement, and the third replacement is the new rightmost representative `r_(k+1)`.

In each case the next untouched original vertex remains adjacent to the rightmost processed vertex, and all processed vertices lie on one loop. No move splits a loop. Since every loop after initial closure contains an original suffix vertex, after stage `w-1` all loops have merged. There are at most `w-1` moves, and at most two vertices are added per move, so the final meander has order `m+w/2+t` for some `0<=t<=w-1`.

**Why the stage/type record determines the inverse.** The record is a word of length `w-1` over `{skip,CO,OO,CC}`. Work backward from the output, undoing later stages first. Immediately before undoing stage `k`, all later moves have been undone, while the insertions from stages `<k` remain. Let `t_<k` be the number of OO/CC letters among those earlier stages, known from the record. Then the processed prefix just before stage `k` has exactly `p=k+2t_<k` vertices. In zero-based suffix coordinates, the marked short upper arch for Lemma 5's inverse is at position `p-1` for CO or CC, and at position `p+1` for OO. The latter follows because the mirror wiggle replaces the old rightmost processed vertex by three, and the CO exchange is at its third vertex. A `skip` requires no inverse. Lemma 5 gives a unique inverse from this marked position and case. After all stages have been undone, delete the first `w` added vertices and restore their incidences as open stubs. Thus the output meander and the four-letter word uniquely determine the input half-state. There are only `4^(w-1)` words, with no independent position marks. A word `a` also fixes `t(a)`, its number of OO and CC letters, so the output order is exactly `m+w/2+t(a)`. Hence the injective map to `(a,output)` gives the stated sharper half-state count without summing over all possible final orders for each word.

**Full cut count.** There are at most `(w+1)^2` admissible nonnegative even pairs `(i,j)` with `2<=i+j<=w`. For a fixed pair the left and right half-state bounds multiply; the noncrossing cross-gap gluing is unique, and restricting to connected full meanders only decreases the count. The `t` sum contributes at most `w` and the size shifts cost at most `C^(3w/2)` in each half. Therefore a crude bound is

`(w+1)^2 * 16^(w-1) * C^(3w) M_m M_n <= (64C^3)^w M_mM_n`,

using `(w+1)^2<=4^w` for `w>=2`. The width `w` need not be fixed as `m,n` grow. The record-location and sharper count were independently checked in `1342_alpha_tree_entropy/REFEREE_LINEAR_RECORD_COMPLETION.md`. ∎

**Bridge corollary.** For an even order `N=2m`, the centered atomic off-block event `E_delta` satisfies `F_(N,K)` with `K<=ceil(20Ndelta+20)`, and `F_(N,K)` implies at most `2+2K` cut arches. The exact junction lies in `E_delta`, so

`q_(N/2) <= P_N(E_delta) <= exp(O(Ndelta+1)) q_(N/2)`, `q_(N/2)=M_m²/M_(2m)`.

Thus for `delta_N=ell_N/N`, the terminal loss is `N^o(1)` whenever `ell_N=o(log N)`, a wider resolution window than the previous `ell_N log(ell_N+1)=o(log N)` condition. This still does not prove the continuum small-ball exponent or rare-event probability transport.

**Exact finite diagnostic.** A direct implementation of the initial closure and this fixed-stage CO/OO/CC schedule checked every one-loop meander at orders 4 and 6 (42 and 1,828 systems). At every stage the pair was adjacent, the resulting matchings remained noncrossing, and the final system was one loop. The inverse-location argument above is analytic and does not depend on this check.
