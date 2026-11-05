# What the certificate verifies (my statement, written before reading the project proofs closely)

Letters: 1 = step +, 0 = step -. A state is (s,h): a word s and a head position h.

## The machine A (Proposition T of 637)

Outgoing moves of (s,h), n = len(s):
* L: (s,h-1) if h>0; R: (s,h+1) if h<n-1.
* split: if s[h]=1, replace it by 101; if s[h]=0, replace it by 100 or by 001; new head h+1 (the middle letter, always 0).
* merge: if 0<h<n-1 and the triple s[h-1:h+2] is 010, 011 or 110 (middle letter 1), replace it by 0, 1, 1; new head h-1.
  The triple 111 (+ [+] +) has NO merge.
Root: (101, 1). `tape_count.py` confirms: closed walks root -> root of length 2m-2 number 1, 2, 8, 42, 262, 1828, 13820,
110954, 933458 for m = 1..9, the meander numbers M_m. (Independent of the proof of Proposition T.)

## The finite SOS statement (T2)

Let X be any finite set of states, A_X the 0/1 adjacency of A restricted to X (principal submatrix), W, Z > 0 vertex weights,
D_i = sqrt(Z_i/W_i), B = D A D^-1 (>= 0), S = (B + B^T)/2. For every unordered pair {i,j} choose f_ij > 0 with f_ji = 1/f_ij.
Row value
    c_i = (1/2) sum_{j: i->j} f_ij W_j/W_i + (1/2) sum_{j: j->i} f_ij Z_j/Z_i          (both sums over j in X).
If c_i <= R for all i in X then lambda_max(S) <= R, so rho(A_X) <= R.
Proof (mine): with t_ij = f_ij u_j/u_i, u = sqrt(WZ), one checks S_ij t_ij = (1/2)(A_ij W_j/W_i + A_ji Z_j/Z_i) f_ij, so c_i = sum_j S_ij t_ij.
t_ij t_ji = 1 gives 2 v_i v_j <= t_ij v_i^2 + t_ji v_j^2, hence v^T S v <= sum_i c_i v_i^2 <= R |v|^2. For a Perron vector v >= 0 of B
(B v = rho v) one has v^T S v = v^T B v = rho |v|^2. [rho(B) = rho(A) since B is similar to A.]

## From rows to mu (normalization)

Let n = 2m-2 and X_n = all tapes (in a strip of height 2m, with length <= n+3). Every closed walk of length n from the root stays in X_n,
so M_m = (A_X^n)_{root,root}. Let C be the strongly connected class of the root and v >= 0 the Perron vector of B restricted to C
(positive on C, zero elsewhere, B v = rho_C v with rho_C <= rho(B) and v^T S v = rho_C |v|^2 <= R|v|^2).
Then (B^n)_{rr} v_r <= (B^n v)_r = rho_C^n v_r, so M_m <= rho_C^(2m-2) <= R^(2m-2), and mu = lim M_m^(1/m) <= R^2.
So R is the bound per STEP of the walk (one step = half a crossing pair), and mu <= R^2 (not R, not R^4).

## From infinite state space to finitely many formal rows (padding)

Weight family (r, k, q, d): for a padded word s' (the actual word with r+1 zeros added at both ends), head h',
    W(s,h) = U[s'[h'-r .. h'+r]] * prod_{k-grams g of L'} F[g] * prod_{k-grams g of R'} F[g],
L' = s'[:h'] read left to right, R' = s'[h'+1:] read right to left;  Z(s,h) = same formula on the complemented padded word.
Head shifts: f_ij = E[cut] (right, cut = s'[h'-q+1 : h'+q+1]) and f = 1/E[s'[h'-q : h'+q]] (left; same cut for the reverse edge).
Edit edges (short tape with head sh, long tape with head sh+1 differing by a triple at the head): f(short -> long) = G[window of radius d
around the SHORT head in the padded short word, long triple], f(long -> short) = 1/G.
The three classes of edges are exactly the support of S (head shifts, splits, merge-reverses); the only pair of states joined by NO edge
that looks like an edge is (+ [+] +, merged letter), which A does not contain (and the transposed split also does not exist).
Claim: for every actual state i of any restriction X, c_i (restricted to X) <= formal row of the window of 2r+3 letters of the padded
word around the head, where the formal row is c computed on that window as a stand-alone word with ALL neighbours. Reasons:
(1) restriction only deletes positive terms; neighbours that need letters beyond the actual tape are present in the formal row;
(2) every ratio W_j/W_i, Z_j/Z_i, f_ij depends only on letters within distance max(r+1 (cut q, grams k), r+1 (edit window d+1)) of the head
    (grams far from the head cancel), and the formal row has r+1 letters on each side; for a merge the shortened word still has the short window.
Hence all c_i <= R follows from the finite check over all 2^(2r+3) words ("formal rows"). The verifier's 404 physical states and 2,765
inverse checks are implementation sanity checks, not part of the proof; my own literal test (spot.py) does the same job with exact arithmetic
on long tapes.

Then: all formal rows <= R => c_i <= R on every finite restriction => rho <= R => M_m <= R^(2m-2) => mu <= R^2.

(Erratum to the remark in the padding section: there is simply no edge for the triple 111; a state whose head triple is 111 has no
merge term, outgoing or incoming, in A, and the formal row omits it too. Nothing is "deleted on one side"; the row of the true graph lacks the term.)
