# Removing the logarithm from the unconditional counting loss

Status: proved, independently reviewed in `1341_alpha_tail/REFEREE_SQRT_COUNT_BOUND.md`, including the ratio-free assembly. No claim of the conjectured power law.

## Theorem

Let M_n count rooted connected meanders of order n, and let mu=lim M_n^(1/n). There is an absolute K such that

    mu^n exp(-K sqrt(n)) <= M_n <= mu^n

for every n>=1. This improves the exp(-K sqrt(n) log n) lower factor of 628. The constant is not evaluated, and no new numerical bound for mu is asserted.

Inputs: junction-product supermultiplicativity from 120; the elementary Catalan-pair bound mu<=16; `PROOF_SINGLE_NOOSE.md`; `PROOF_DISC_COMPLETION.md`; the shallow-cap lemma in 1343. The geometric lemmas use ordinary planar topology and the named sphere-cut decomposition theorem. Every counting factor below is retained before being absorbed. The main theorem does not require the ratio theorem of 113.

## Balanced reverse inequality

For sufficiently large N choose the single balanced separator. Put W=ceil(4 sqrt(27N)). Its number w of ports satisfies w<=W. The two disc patches have original crossing counts v_1+v_2=2N, with N/2<=v_i<=3N/2. Complete them by the recognizable-cycle construction. Let j_i be their completed orders. They obey

    v_i/2 <= j_i <= v_i/2+4w.

The two rooted completed meanders, their protective-cycle turn words and interior-side bits recover both patches. Their cyclic port orders glue them back together: an orientation and a cyclic alignment suffice, at most 2w choices. Invalid port color matchings can simply be discarded. Root and orientation of the original meander can be restored with at most O(N) choices. The combinatorial embedding of a disc patch and cyclic boundary identification determine the glued colored spherical map; there is no metric or continuous gluing datum. Homeomorphisms of a circle with the same port alignment introduce no additional combinatorial possibilities.

A deterministic choice of one valid separator and completion for each original meander therefore gives an injection into these records. Equivalently, the described reconstruction maps all valid records onto the original meanders, which is enough for an upper bound. The record count is at most a polynomial in N times A^(2W) times the largest product M_(j_1)M_(j_2). Even if both completed meanders retain arbitrary original root conventions, marking one ring vertex in each costs at most O(N)^2 in total. The polynomial also accounts for w, the two output sizes, cycle lengths if counted separately, alignment, and the original root. All output sizes are O(N), so these are only a fixed number of polynomial factors; no exponent of N depends on w.

Consequently there is a constant B with

    M_N <= exp(B sqrt(N)) max M_a M_b,                       (1)

where the maximum ranges over integer pairs satisfying

    N<=a+b<=N+8W,   N/4<=a,b<=3N/4+4W.

The polynomial record factors have logarithm O(log N), hence are absorbed in B sqrt(N). The original crossing counts v_i may be odd; only the completed orders a,b must be integers. For all sufficiently large N the displayed size bounds imply a,b in [theta N,(1-theta)N], with theta=1/10.

## From balanced splits to all sizes

Define L(n)=n log(mu)-log(M_n). Supermultiplicativity and Fekete's lemma give L(n)>=0, L(n)=o(n), and subadditivity L(m+n)<=L(m)+L(n). Equation (1) gives, for each sufficiently large N, some a,b in [theta N,(1-theta)N], theta=1/10, such that

    L(a)+L(b) <= L(N)+(a+b-N)log(mu)+B sqrt(N)
               <= L(N)+B_1 sqrt(N).

Here the added size is O(sqrt(N)) and mu<=16. Taking the weighted average with denominator a+b>=N, and using L(N)>=0, at least one child k in {a,b} satisfies

    L(k)/k <= L(N)/N+B_1/sqrt(N).                            (2)

Fix a large target scale n and start from an arbitrarily large integer N_0. Repeatedly choose a child satisfying (2), until the first size k<=n. Then theta n<k<=n, and every step shrinks by a factor between theta and 1-theta. The accumulated sum of reciprocal square roots of the ancestor sizes is at most C_theta/sqrt(n): going backward from the last ancestor, sizes grow at least geometrically with ratio 1/(1-theta). Thus

    L(k)/k <= L(N_0)/N_0 + B_1 C_theta/sqrt(n).

Choose N_0 sufficiently large that its first term is at most 1/sqrt(n). It follows that for every sufficiently large n there exists k in [theta n,n] with L(k)<=B' sqrt(n), uniformly in n.

For any target n, remove such a block k and repeat on the remainder r=n-k. Each remainder is at most (1-theta) times its predecessor. Subadditivity bounds L(n) by the sum of the block losses plus a bounded terminal remainder. Their square-root bounds sum as a geometric series, giving L(n)<=K sqrt(n). Enlarge K to cover the finitely many small n. Exponentiating proves the theorem. The upper bound M_n<=mu^n already follows from supermultiplicativity.

## Scope

This is an unconditional finite-size counting improvement; it needs no ratio theorem, scaling limit, SLE estimate, field nonconcentration, or probability transfer. It does not establish a power-law prefactor, the value of alpha, or positivity. The protected cycle is essential: without retaining it, the cap construction alone leaves the old N^O(w) location ambiguity.
