# 641: proofs for the tape-window bounds

Status words as in `WAVE_620_COMMON.md`. Everything in this file is conditional on **Proposition T** of
`637_rung2_matrix_product/PROOF.md` section 1 (proved there, not refereed), which rests on Theorem M of
`630_referee_band_growth/REPAIRS.md` (refereed). Nothing else from the project is used. M_m is the number of
meanders of order m and μ = lim M_m^(1/m).

## 0. The machine (quoted from 637)

A *tape* is a finite word over {+, −} with a *head* on one letter; I write … x [s] y …. Moves:

    left    head to the previous letter (if there is one)        right   head to the next letter (if there is one)
    split   [+] → + [−] +        [−] → + [−] −        [−] → − [−] +
    merge   − [+] − → [−]        − [+] + → [+]        + [+] − → [+]      (the head needs both neighbors; no merge at + [+] +)

The *root* is + [−] +. **Proposition T (637, no cap):** M_m is the number of walks of length 2m − 2 from the root
to the root.

**Class C.** A tape is of class C if its first and its last letter are + and it has one more + than −.
*The class is closed under the moves and contains the root* [proved]: left and right do not change the word; a
split replaces s by three letters with the same sum, and if s is the first letter then s = + (the only split is
+ − +) or s = − is not the first letter; in all cases the first and last letters stay +. A merge replaces x + y by one
letter with the same sum; if x is the first letter then x = +, so y = − and the new letter is +; the same at the
other end. So every tape reached from the root is of class C. (`machine.py` confirms that the tapes reached from the
root with at most J pairs are exactly the class-C tapes, J ≤ 8; this is not used.)

A tape with 2j + 1 letters has *j pairs*. A walk root → root has even length (the number of letters left of the
head changes by ±1 in every move).

## 1. The two comparison lemmas [proved]

**Lemma U (upper).** Let W > 0 on the tapes of class C and R > 0 with Σ_(T → T′) W(T′) ≤ R·W(T) for every tape T of class
C (sum over all moves). Then M_m ≤ R^(2m−2) for all m ≥ 1, so μ ≤ R².

*Proof.* Let c_n(T) be the number of walks of length n from the root to T. By induction
Σ_T c_n(T) W(T) ≤ R^n W(root): Σ_T′ c_(n+1)(T′) W(T′) = Σ_T c_n(T) Σ_(T → T′) W(T′) ≤ R Σ_T c_n(T) W(T). All terms are
non-negative, so c_n(root) W(root) ≤ R^n W(root). Proposition T. ∎

**Lemma R (return).** From a tape of class C with j pairs the root is reached in at most j(2j + 1) + 2 moves, and
no tape on the way has more than j pairs.

*Proof.* Let the tape have N = 2j + 1 ≥ 5 letters. The N − 2 inner letters are j letters − and j − 1 ≥ 1 letters +,
so two neighboring inner letters differ; let p be the position of the + of such a pair. It has both neighbors and
one of them is −, so after moving the head to p (at most N − 1 moves, each allowed, the word does not change) the
merge is allowed. The result has N − 2 letters and is of class C. For N = 3 the word is + − + and at most one head
move gives the root. ∎

**Lemma L (lower).** Let u ≥ 0 on the tapes of class C, zero on every tape with more than J pairs, u(root) > 0, and
R > 0 with Σ_(T → T′) u(T′) ≥ R·u(T) for every tape T of class C with u(T) > 0. Then μ ≥ R².

*Proof.* Let S = {u > 0} and let w_n(T) be the number of walks of length n from the root to T that stay in S. Then
Σ_T w_n(T) u(T) ≥ R^n u(root) by induction (successors outside S contribute 0 to the left side of the hypothesis).
S is finite, so w_n := Σ_T w_n(T) ≥ R^n u(root)/max u. Fix for every tape of S one return path of Lemma R, of length at
most D = J(2J + 1) + 2. Appending it is injective on the walks with a given end tape and a given total length, so
Σ_(d = 0..D) c_(n+d)(root) ≥ w_n. Hence max over d ≤ D of M_((n+d)/2 + 1) ≥ R^n u(root)/((D + 1) max u), and letting
n → ∞ gives μ ≥ R². ∎

Lemma R is the statement that 606 had as "data only" for the label machine; in tape form it is two lines.

## 2. The weight family and the row lemma [proved]

Fix integers r ≥ k ≥ 1, positive numbers F[c] for c ∈ {0, …, 2^k − 1}, and positive numbers U[w] for
w ∈ {0, …, 2^(2r+1) − 1}. Letters are coded + = 1, − = 0; a word of letters is read as a binary number with its
first listed letter as the most significant digit.

- *k-grams of a tape.* For every letter left of the head that has at least k − 1 letters before it: the k letters
  ending at it, listed left to right. For every letter right of the head that has at least k − 1 letters after it:
  the k letters starting at it, listed right to left. (Both lists are "far to near"; the family is mirror symmetric.)
- G(T) = product of F[c] over the k-grams c of T.
- The *window* of a tape whose head has at least r letters on each side is the word of the 2r + 1 letters from r
  places left of the head to r places right; V(T) = G(T)·U[window of T].
- A *row* is a word ω of 2r + 3 letters, regarded as a tape with the head on its middle letter.
- **Row inequality with constant R:**  Σ over the moves ω → ω′ of (G(ω′)/G(ω))·U[window of ω′] ≤ R·U[window of ω].

**Lemma 1 (locality).** Let T be a tape whose head has at least r + 1 letters on each side, T = α ω β with ω the row
around the head. The moves of T are T → α ω′ β for the moves ω → ω′ of ω; the head of ω′ has at least r letters on
each side; the window of α ω′ β is the window of ω′; and G(α ω′ β)/G(T) = G(ω′)/G(ω).

*Proof.* The head of ω has both neighbors, and a move is determined by the head letter and its two neighbors and
changes only these three letters and the head position by at most one. So the moves of T and of ω correspond, and
the letters of ω at distance at least 2 from its head, in particular its first r and its last r letters, are
unchanged and stay on their side of the head. After a move the head has at least r letters of ω′ on each side
(left or right move: r and r + 2; split: r + 2 and r + 2; merge: r and r), so the window lies in ω′.
k-grams: a letter of α, or of β, has its k-gram inside α, or β, in T and in α ω′ β alike. A letter of ω (or ω′) left
of the head at place i ≥ k of the row (counted from 1) has its k-gram inside the row: these are exactly the k-grams
of G(ω) (or G(ω′)) on the left. A letter at place i ≤ k − 1 has a k-gram made of letters of α and the first i
letters of the row; since k − 1 ≤ r − 1, these letters are the same in ω and ω′, they are left of the head in both,
and they carry no k-gram in G(ω) or G(ω′). The right side is the mirror image. So G(T) = G(ω)·X and
G(α ω′ β) = G(ω′)·X with the same X > 0. ∎

**Lemma 2 (padding).** For a tape T of class C let pad(T) be T with r + 1 letters − put in front and r + 1 letters −
put behind (head on the same letter). Every move T → T′ is a move pad(T) → pad(T′).

*Proof.* Left, right and split: the same move of pad(T). Merge: allowed in T iff the head letter is +, has both
neighbors in T, and they are not both +; then the same holds in pad(T), with the same result. ∎

**Theorem U (conditional on Proposition T).** If every row satisfies the row inequality with the constant R, then
M_m ≤ R^(2m−2) for all m ≥ 1, and μ ≤ R².

*Proof.* W(T) := V(pad(T)) is positive on class C. Let ω be the row of pad(T) (its head has at least r + 1 letters on
each side). By Lemma 2, Lemma 1 and the row inequality,

    Σ_(T → T′) W(T′) ≤ Σ over the moves pad(T) → S of V(S) = G(pad(T))/G(ω) · Σ_(ω → ω′) G(ω′) U[window of ω′]
                    ≤ R · G(pad(T)) · U[window of ω] = R·W(T).

Lemma U. ∎

Remarks.
- There is no worst case anywhere: a row determines the window of every successor and the exact ratio of the
  global factor. (In the bracket-word method the far flip of a CC move is charged its worst case.)
- Every row occurs in a tape of class C (put letters around ω so that the first and last are + and the sum is 1),
  so the reduction to rows loses nothing on class C: the smallest R of the row system is the smallest R for
  which V satisfies the inequality on all tapes of class C whose head has r + 1 letters on each side. (Whether
  every tape of class C is reached from the root is not needed; it holds for at most 8 pairs by `machine.py`.)
- Explicit ratios (used by `tapeup2.c` and `cert1.c`). With LK = code of l_k … l_1 (the k letters left of the head,
  nearest last), RK = code of ρ_k … ρ_1 (the k letters right of the head, nearest last), s the head letter:
  left removes the k-gram LK and adds (ρ_(k−1) … ρ_1 s); right removes RK and adds (l_(k−1) … l_1 s); a split
  [s] → x [−] y adds (l_(k−1) … l_1 x) and (ρ_(k−1) … ρ_1 y); a merge removes LK and RK. `verify_tape.c` and
  `verify_tape.py` do not use these formulas: they list the k-grams of ω and ω′.

## 3. What a certificate is, and what the verifier checks [certified computation + Theorem U]

A certificate is (r, k, B, P, Q, N[0 … 2^k − 1], T[0 … 2^(2r+1) − 1]) with positive integers N[c], T[w]. It claims the row
inequalities with F[c] = N[c]/2^B, U[w] = T[w], R = P/Q. `verify_tape.c` enumerates all 2^(2r+3) rows as words, applies
the moves of section 0 to the word, reads the windows off the letters, lists the k-grams before and after, and
checks Σ (ratio)·T[window′] ≤ (P/Q)·T[window] in unsigned 128-bit integers after multiplying by the least common
denominator; every product and sum is guarded against overflow (a failed guard rejects the certificate).
If it prints PASS, Theorem U gives μ ≤ (P/Q)², conditional on Proposition T.

## 4. Lower-bound tables with the side lengths as index [row reduction: proved in the direction needed; values: data]

Weight u(T) = G(T)·U(window with boundary marks, a, b), a = letters left of the head, b = letters right, zero when
a + b > 2J; left view = last r letters of "^" + (left word), right view = first r letters of (right word) + "$"
(a letter with fewer than k − 1 letters beyond it carries no k-gram; when a boundary is visible this is decided by
the view). A row is the window with one more symbol revealed on each side that does not yet show its boundary.
As in Lemma 1 a row fixes the views of every successor, the exact ratio of G, and the new side lengths:
left (a − 1, b + 1), right (a + 1, b − 1), split (a + 1, b + 1), merge (a − 1, b − 1).

**Admission rule.** The row inequality is imposed at (row, a, b) iff: a view that shows its boundary has exactly
a (or b) letters and its outermost letter is + (a head on the first or last letter is +); a view that does not
show its boundary has nh = a − (its length) ≥ 0 hidden letters, and if nh = 0 its outermost letter is +; a + b is
even; a + b ≤ 2J; and 1 − (sum of the visible letters, head included) lies between loL + loR and hiL + hiR, where
a side with nh ≥ 1 hidden letters has lo = −nh + 2, hi = nh and a side with none has lo = hi = 0.

*Every (row, a, b) realized by a tape of class C is admitted* [proved]: the side lengths and the boundary
conditions are read off the tape; the hidden letters of a side with nh ≥ 1 start (or end) with a + and so sum to a
number between −nh + 2 and nh; the total sum is 1. This is the direction Lemma L needs (an admitted row that is
not realized only weakens the bound). The converse also holds (choose the hidden letters; sums of that parity in
the range are all attained), so no superfluous row is imposed; it is not needed.

With this rule a table whose every admitted row inequality holds with constant R, and whose entry at the root is
positive, gives μ ≥ R² by Lemma L. No such table was certified in this folder (float values only, `lower.py`).
