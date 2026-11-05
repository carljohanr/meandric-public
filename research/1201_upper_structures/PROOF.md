# Gram length up to r+1 in the paired reciprocal certificate

Status: proved (written argument below), not independently reviewed. It extends the weight family used in
`../1007_edit_edge_1230/PROOF.md` by allowing the gram length k to be r+1 instead of at most r. Nothing else in
the family changes. The theorem it rests on is the finite SOS theorem of `../1003_paired_path_consistency/PROOF.md`
section 1 together with the head and edit factor assignments of `../1007_edit_edge_1230/PROOF.md`. Both are used
unchanged and are not re-proved here.

## Notation

A headed tape is (s, h): a word s over {0,1} and a head position h. The left side is L = s[0..h-1], the right side
is R = reverse(s[h+1..]); both are read from the far end toward the head. The weight is

    W(s,h) = U(s[h-r..h+r]) * prod_{k-grams g of L} F(g) * prod_{k-grams g of R} F(g),

evaluated on the tape padded with r+1 zeros on each side. The formal row of a state is the word of length 2r+3 centered
at the head (r+1 letters on each side). The certificate checks every formal row, computing weights on the formal row
word itself (no padding).

## Lemma (locality of the gram part, k <= r+1)

Let k <= r+1. For each of the four forward moves of the tape machine (head left, head right, split, merge) the ratio
of the gram parts of W(t) and W(s) is a function of the letters of s at distance at most r+1 from the head only.
The ratio is the same whether it is computed on the true tape or on the formal row word. The same holds for the
adjoint weight Z(s) = W(complement of s), which is how the incoming terms are evaluated.

Proof. Write the old state (s,h) and the new state (t,h').

* Head left (h' = h-1): L' is L without its last letter s[h-1]; R' is R with s[h] appended at the head end. All
  k-grams of L and R persist except the last k-gram of L, which is s[h-k..h-1]. The k-gram that is new in R'
  ends at the new letter and is (s[h+k-1], ..., s[h]) read toward the head. The letters used are
  s[h-k..h+k-1], at distance at most k <= r+1 on the left and k-1 on the right.
* Head right is the mirror image.
* Split (s[h] replaced by a triple p0 p1 p2, head at the middle, h' = h+1): L' = L followed by p0, R' = R followed by
  p2. New grams end at p0 or p2 and use s[h-k+1..h-1] or s[h+1..h+k-1]. Distance at most k-1.
* Merge (triple s[h-1] s[h] s[h+1] replaced by one letter, head on it, h' = h-1): L' is L without its last letter
  s[h-1] (lost gram s[h-k..h-1], distance k); R' is R without its last letter s[h+1] (lost gram
  (s[h+k], ..., s[h+1]), distance k). The new head letter does not belong to a gram.

In every case grams far from the head are unchanged, so they cancel in the ratio. The data needed lie within
distance k <= r+1 of the head, which is exactly the formal row. The window part U has radius r and the shifted
window radius is at most r around a head that moved by one, so it also lies in the row, as in 1007. For the adjoint
weight the same statement holds after complementing every letter, because Z(s) = W(P s) with P the complement
and incoming moves are the complements of outgoing moves. QED.

## Consequence

All hypotheses of the finite SOS theorem hold for k up to r+1 with the same proof as in 1003 and 1007: the
weights are positive, every actual inverse pair of edges has reciprocal factors, and the row sum of any actual
physical state with boundary or cap restrictions is a sub-sum of the formal row sum of the padded state. The
exact checker `../1007_edit_edge_1230/verify.py` evaluates formal rows on the formal words, which by the lemma are
the true rows, and tests the actual physical states through four pairs against their padded formal rows.

The earlier code limited k <= r because its evaluator class rejected larger k; the theory never needed it.

## Small-radius check (data)

At a = r+1 = 5 the converged family optimum with k = 5 is 12.36392 against 12.36583 for k = 4 (log-barrier
Newton solver `ipm.py`), at a = 6 it is 12.33556 against 12.33585. These are numerical optima; the certified
values of any new witness come from `verify.py` only.
