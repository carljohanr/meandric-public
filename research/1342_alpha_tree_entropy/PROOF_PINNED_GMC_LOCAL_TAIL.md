# Pinned GMC inner-mass band and annulus localization

**Status.** Proved for a zero-boundary GFF in a fixed planar disk, and for a residual field obtained by subtracting any independent smooth rank-one Gaussian scalar. This is a local input, not an identification of the unit-area quantum-sphere two-clump probability.

Let `D` be a disk, `z0` an interior point, `gamma in (0,2)`, and `h` its zero-boundary GFF with covariance `G_D(x,y)=-log|x-y|+f_D(x,y)` in the chosen normalization. On a compact disk `U` about `z0`, `f_D` is smooth. Suppose `h=B phi+R`, where `B` is a nondegenerate centered Gaussian scalar independent of `R`, and `phi` is smooth and equals one in a neighborhood of `z0`. Such a decomposition can be obtained by projecting onto a smooth Cameron–Martin direction and rescaling it; the covariance of `R` is `-log|x-y|+f_R(x,y)` with `f_R=f_D-Var(B)phi(x)phi(y)`, still smooth on `U×U`. Write `M_gamma^R` for the usual subcritical GMC. Put `Q=2/gamma+gamma/2`, `q=2(Q-gamma)/gamma=4/gamma²-1`, and `a=4/gamma²=q+1`.

For a sufficiently small fixed `r>0`, define the pinned inner mass

`A_r = int_{B_r(z0)} |x-z0|^{-gamma²} g(x) M_gamma^R(dx)`,                 (1)

where `g` is a strictly positive smooth function on `U`; this includes the deterministic smooth factors from a gamma insertion. Then there is `c_r in (0,infinity)` such that

`P(A_r>t) ~ c_r t^{-q}` and `P(t<=A_r<=2t) ~ c_r(1-2^{-q})t^{-q}`.       (2)

The same conclusion holds after multiplying `A_r` by any positive random scalar independent of `R` with a finite `(q+eta)` moment (in particular `exp(gamma B)`), by Breiman's lemma; if `B` is conditioned to a fixed compact interval the two-sided band estimate remains of order `t^{-q}` directly.

**Primary theorem and tail transfer.** Wong, *Universal tail profile of Gaussian multiplicative chaos*, arXiv:1902.04054, Proposition A.1 proves the exact-kernel inserted mass `int |x-z0|^{-gamma alpha}M_gamma(dx)` has positive regularly varying tail with exponent `2(Q-alpha)/gamma` for `alpha in (gamma/2,Q)`. Here `alpha=gamma` is admissible because `gamma<2`, and gives `q`. Proposition A.1 alone applies to the exact reference kernel; the smooth-covariance transfer uses Wong's Lemma 3.8, equation (3.16), whose interpolation is expressly for `-log|x-y|+f(x,y)` with the decomposition hypothesis. A smooth `f_R` on a compact subdisk satisfies that hypothesis by Wong Lemma 2.3. Wong defines the pinned functional `M_{gamma,g_W}(z0,A)=int_A |x-z0|^{-gamma²} exp(gamma² f_R(x,z0)) g_W(x) M_gamma^R(dx)` in Lemma 3.1. Choose `g_W(x)=g(x)exp(-gamma² f_R(x,z0))`; then this functional is exactly (1), and `g_W(z0)>0`. Lemma 3.8, together with the localization Lemma 3.7, gives a finite positive limit for `lambda^a E[A_r^{-1} exp(-lambda/A_r)]`. Apply Karamata's theorem to the measure `v P(A_r^{-1} in dv)` near zero and then integrate `v^{-1}`: this is equivalent to `P(A_r>t)~c_r t^{-(a-1)}` with `c_r>0`. Equivalently, this is the pinned counterpart of Wong's Corollary 2.11. The fixed-factor band formula is subtraction of the tails at `t` and `2t`; no density or local limit theorem is needed.

**Annulus.** Fix `0<r0<r1<dist(z0,partial D)` and let `V={r0<|x-z0|<r1}`. Let `A_V=int_V w(x) M_gamma^R(dx)` where `w` is a bounded positive smooth weight; this includes the gamma-insertion weight `|x-z0|^{-gamma²}g(x)` because it is bounded on `V`. Wong Theorem 1.1 / Corollary 1.2 applies to `R` on a compact neighborhood of `V`, so `P(A_V>t)=O(t^{-a})=o(t^{-q})`. Thus for each `epsilon>0`,

`P(t<=A_r<=2t, A_V<=epsilon t) = c_r(1-2^{-q})t^{-q}+o(t^{-q})`.  (3)

No independence between inner mass and annulus mass is used: subtract `P(A_V>epsilon t)=O(t^{-a})` from the band probability.

**Remaining sphere step.** Equations (2)–(3) do not give the probability of two balanced clumps in a *unit-area quantum sphere*. That requires an explicit law or absolute-continuity comparison for the sphere's fixed conformal embedding and global area normalization, including any zero-mode/Palm weighting and the second clump. The local tail is robust to smooth finite-rank Gaussian subtraction, but this alone does not transfer a `t^{-q}` event to the normalized sphere.

Primary source: https://arxiv.org/pdf/1902.04054, Proposition A.1, Lemmas 2.3, 3.7, 3.8, Theorem 1.1 and Corollary 1.2.
