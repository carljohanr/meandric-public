# 808: repairs to 805 PROOF.md

None of these changes a statement, a constant or the certificate. Items 1 and 2 are rewordings inside the proof of
Theorem F; item 3 completes a remark that is not used; item 4 corrects an equality clause in section 4.

## 1. Step 1, the sentence before the case list

Printed: "In each case the letters of ω′ outside W′, given (V, m), are chain continuations from the two outer
letters of W′ on each side."

This is not literally true for m = L (and R): there W′ = f2 f1 w, and the letters g1, g2 right of W′ are part of V,
so given (V, m) they are not random.

Replace by: "In each case the conditional law of the letters of ω′ outside W′, given (W′, m), is the kernel
K(W′, ·) = (the chain continued to the left from the first two letters of W′) ⊗ (the chain continued to the right
from the last two letters of W′). Indeed, under ν and given the letters of ω at any set of consecutive places
containing the window, the letters outside are continued by the chain from the two end pairs (order-2 Markov
property, two sides independent), and the label is drawn from p(· | w), a function of the window alone; and in
each of the three cases below W′ is a function of (m and) a set of consecutive letters of ω containing the
window, the letters of ω′ outside W′ are the letters of ω outside that set in the same order, and the end pairs of
the set are the end pairs of W′." After the case list add: "K does not depend on m. Hence the law of ω′ is
(law of W′) ⊗ K, where the law of W′ is the left side of (S)."

Check of the three cases: m = L: the set is f2 f1 w (places −r−2, …, r), W′ = the same word. m a split: the set is
the window w, W′ = w_1 … w_r x′ − y′ w_(r+2) … w_n, whose end pairs are w_1 w_2 and w_(n−1) w_n because r ≥ 2. m = M: the
set is V, W′ = f2 f1 w_1 … w_(r−1) z w_(r+3) … w_n g1 g2, whose end pairs are f2 f1 and g1 g2.

## 2. Step 4, the sentence on the finite tape

Printed: "Along each of these sequences the head keeps at least two letters of T_0 on each side, and moves read
and change only the head letter and its two neighbors, so each sequence is a walk of the tape machine from T_0."

The letters of T_0 themselves may be merged away; what is kept is a number of letters.

Replace by: "Run the labels on the finite tape T_0. By induction on i, after i moves the finite tape is the part
of ω_i from ℓ + A_i letters left of the head to ℓ + B_i letters right of it, and on E both numbers are at least
ℓ − εn ≥ 2: a move reads and changes only the head letter and its two neighbors, which lie inside this part, and
it changes the two numbers by (Δa, Δb). In particular every label applicable at ω_i is applicable at the finite
tape (L needs a letter on the left, R one on the right, M both neighbors). So each sequence is a walk of the tape
machine from T_0."

Also, for the record, the count in the same step: the build walk has at most 3(4ℓ + 5) ≤ 12εn + 51 moves and the
return walk at most 3(4ℓ + 5 + εn) ≤ 15εn + 51, so n′ ≤ 27εn + 102, inside the printed 30εn + 110.

## 3. Remark 3 of section 1 (strips): the three adaptations

For a strip of K levels the state is (tape, level of one vertex), the far tape is a stationary nearest-neighbor
Markov path in the strip, and a split is allowed iff the new path stays in the strip. Changes to the proof:
(i) *Lemma G.* The path chain is irreducible with period 2 and has forced steps at the walls, so "positive 3-gram
probabilities" fails. Irreducibility is enough: the stationary measure of an irreducible finite chain is
shift-ergodic. The density of ν with respect to it is positive on all windows that lie in the strip, by (P).
(ii) *Lemmas R, B.* A merge and a backward split remove two vertices of the path and leave the others at their
levels, so both sweeps stay in the strip; for B the forward walk passes only through tapes whose vertex levels
are among those of the target.
(iii) *Padding.* T_0 comes with its levels. Put in front the letter + (from level −1 to 0 in the coordinates of 637)
followed by the monotone path from 0 to the first level of T_0, and behind it the monotone path from its last
level to −1 followed by +. Both stay in the strip, and at most 2K + 2 letters are added.
With these the proof of Theorem F goes through [sketch; I did not write out more than this]. The values of the
remark are confirmed by my own program (REVIEW.md, section 3).

## 4. Theorem D (b), the equality clause

Printed: "with equality iff φ∘c′ is constant on the moves of the binding rows."

For each row ρ, Jensen gives log R′_(Ue^φ)(ρ) ≥ log R′_U(ρ) + Σ_m p^U_ρ(m) φ(c′(ρ, m)) − φ(c(ρ)), with equality iff φ∘c′(ρ, ·) is
constant on the moves of ρ that have positive probability. For finitely many rows, the two infima are equal iff
this holds for at least one row at which the right side attains its minimum (if ρ* is such a row then
min(left) ≤ left(ρ*) = right(ρ*) = min(right) ≤ min(left); conversely a row minimizing the left side then also
minimizes the right side and has equality in Jensen). Replace "of the binding rows" by "of at least one row that
attains the minimum in avg(U, φ) (finitely many rows)". Part (c) uses only the inequality and is not affected.
