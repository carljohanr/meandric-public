# Exact normalization of the APPS determinant tilt

Status: complete argument, independent review pending. The statement concerns the determinant defined through Polyakov–Alvarez in APPS, not an independently established spectral zeta continuation for their nonsmooth metric. Primary source: Ang–Park–Pfeffer–Sheffield, https://arxiv.org/html/2005.11845 , Definitions 6.2, 6.4, 6.7, Lemmas 6.3, 6.5, Proposition 6.8 and Theorem 6.9.

## 1. The polynomial factor

Fix the deterministic smooth sphere metric g0 from Definition 6.7, flat on a neighborhood of the unit square. Write K0 for its Gaussian curvature and h for the GFF with curvature pairing zero. For a finite terminal partition P with n leaves, put u=h^P/Q and V=Vol_(exp(2u)g0)(sphere). The regularized metric is g=exp(2(u+C))g0 with

    C=(1/4)log n−(3/4)log V.

Let A0=det'_zeta(Delta_g0)/Vol_g0(sphere), a positive deterministic constant independent of P,n,Q. The Polyakov–Alvarez definition, equation (6.7), gives

    log det'_PA(Delta_g)
      = −(1/(12π)) ∫|∇u|² − (1/(6π))∫K0(u+C)
        +2C+log V+log A0.

The curvature pairing of u is zero. Gauss–Bonnet gives ∫K0=4π. Hence the terms outside the energy are

    −2C/3+2C+log V+log A0
      = (1/3)log n+log A0.

With the Dirichlet convention (f,f)_∇=(1/(2π))∫|∇f|², Lemma 6.5 gives ∫|∇u|²=2π Σ_(S∈P) x_S² where x_S=alpha_S/Q. Therefore

    det'_PA(Delta_g)=A0 n^(1/3) D_P(x),
    D_P(x)=exp(−Σ_S x_S²/6).                         (1)

This derives the n dependence rather than discarding it under conditioning. In particular it does not use the displayed equation (6.9), which contains ∇(2h^epsilon/Q). That display is inconsistent by a factor of four with (6.7), Definition 6.7, and the Gaussian tilt in the proof of Theorem 6.9. Formula (1) follows directly from the definition and (6.7).

## 2. Rank and unconditioned densities

The finite projection has dimension n, not n−1. To see this, the test functions before inversion are 1_S−a_S K0, with a_S chosen to make the integral zero. A linear combination vanishing as a distribution restricts on the square interiors to b_S=0 for each S, since K0 vanishes throughout a neighborhood of the unit square. Thus these n functions are independent. The inverse Laplacian on mean-zero functions is injective, and so is their orthonormalization.

For each fixed finite dyadic partition P, its n projection coordinates alpha_S are independent standard normals. They need not be independent of the event P^epsilon=P. Instead, Lemma 6.3 says that this event is a measurable constraint on these coordinates alone; in coordinates x=alpha/Q it is a Q-independent set Lambda_P. Ancestor averages are weighted sums of leaf averages, so no extra infinite-dimensional condition is hidden in the event.

There are finitely many ordered full 4-ary subdivision trees with any fixed number n of leaves, and each determines a spatial partition. Consequently, if p_n(Q)=P_Q(N=n),

    I_n(s)=Σ_(P: #P=n) ∫_(Lambda_P) exp(−s|x|²/2) dx,
    p_n(Q)=(Q/sqrt(2π))^n I_n(Q²).                  (2)

The sets are empty unless n≡1 (mod 3), since each split adds three leaves. Equations here include zero values on the other integers. For Q≥2, finiteness of N almost surely is the APPS subdivision input; if Q<2, the same formulas describe the finite-N part only.

## 3. Exact tilted partition functions

Let k be real with Q_k²=Q²+k/3>0, and define the two unconditioned weights

    Z_n^E(k)=E_Q[D^k;N=n],
    Z_n^PA(k)=E_Q[(det'_PA Delta_g)^k;N=n].

Multiplication of the Gaussian density by exp(−k|x|²/6) changes its precision to Q_k². All integrals are finite under the stated inequality; no exchange of signed sums is needed. Equations (1) and (2) give

    Z_n^E(k)=(Q/Q_k)^n p_n(Q_k),
    Z_n^PA(k)=A0^k n^(k/3)(Q/Q_k)^n p_n(Q_k).    (3)

Thus, if Q_k≥2,

    Σ_(m≥n) (Q_k/Q)^m Z_m^E(k)=P_(Q_k)(N≥n),
    Σ_(m≥n) A0^(−k)m^(−k/3)(Q_k/Q)^m Z_m^PA(k)
      =P_(Q_k)(N≥n).                              (4)

These are equalities of nonnegative sums. For Q_k<2 their right sides must instead restrict to finite N. This proves the full normalization step of 159 for the energy weight, and identifies the missing polynomial factor for the actual PA determinant.

## 4. What would follow from the tail lemma

If p_(≥n)(Q_k)=n^(−x_k+o(1)) with x_k=4/gamma_k² and Q_k=2/gamma_k+gamma_k/2, then (4) has exponent x_k. At base Q²=25/6 and k≥−1/2, algebra gives x_k+1/x_k=(13+2k)/6. This remains conditional on the cell-count tail and gives no local asymptotic or derivative asymptotic by itself.

There is also an exact interpretation of what happens if only the exponential counterweight is used on the raw determinant. Suppose a nonnegative integer random variable has tail t^(−x+o(1)), x>0. For every real a<x,

    E[N^a;N≥n]=n^(a−x+o(1)).                      (5)

For the upper bound when a≥0, sum over dyadic bands and use tail≤t^(−x+delta) with delta<x−a. For a<0, N^a≤n^a on the tail. For the lower bound, restrict to n≤N≤n^(1+eta). Tail asymptotics make the removed tail negligible compared with P(N≥n), choosing delta small relative to eta. On that band N^a≥n^a for a≥0 and ≥n^((1+eta)a) for a<0. Let eta and delta decrease to zero. This proves (5).

Taking a=k/3 therefore yields exponent x_k−k/3, not x_k, for Σ_(m≥n)(Q_k/Q)^m Z_m^PA(k), up to the constant A0^k. For k≥−1/2, a<x_k follows from k/3=x_k+1/x_k−13/6 and x_k≥1. The difference is a normalization issue, not a contradiction to the fixed-n central-charge tilt.

## 5. Remaining gaps

The actual APPS tail exponent is not proved here. A tail at fixed parameter also does not justify differentiating an o(log n) error in that parameter, so the proposed tree-entropy coefficient needs separate control. The local law must be formulated on the attainable lattice n≡1 (mod 3), and cannot follow from a tail estimate alone without a regularity argument.
