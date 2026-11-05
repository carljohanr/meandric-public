# Sharp logarithmic upper bound for the APPS cell-count tail

Fix ε>0 and Q>2. Let h be the curvature-normalized sphere GFF used by Ang–Park–Pfeffer–Sheffield (APPS), Definition 6.2, and N the number of terminal dyadic squares. Put

    γ=Q−√(Q²−4),     x=4/γ².

**Theorem.** For every real 1≤p<x, E[N^p]<∞. Consequently

    limsup_{n→∞} log P(N≥n)/log n ≤ −x.

This proves the upper inequality of the proposed logarithmic tail at every c<1. At c=−4, x=(17+√145)/12. It does not prove the matching lower inequality.

**Proof.** Write L=log 2. Let I_j count all dyadic squares of depth j satisfying the APPS active condition e^{h_S/Q}2⁻ʲ>ε, even if an ancestor was inactive. The number of split squares in the actual subdivision is at most I=Σ_{j≥0}I_j, and N≤1+3I. For any a>0, 1_{A>ε}≤(A/ε)^{aQ}. Therefore

    I_j≤ε^(−aQ)e^(−aQjL) Σ_{depth S=j}e^{a h_S}.

The square-average covariance of the normalized sphere GFF has the local form −log|z−w|+H(z,w), with H bounded and smooth on a neighborhood of the unit square. In particular Var(h_S)=jL+O(1) uniformly. Define the normalized chaos approximation

    M_{a,j}=4⁻ʲΣ_{depth S=j} exp(a h_S−a²Var(h_S)/2).

Then a constant C depending on a,Q,ε gives

    I_j≤C exp((2−aQ+a²/2)jL) M_{a,j}.            (1)

We justify the standard uniform moment estimate

    sup_j E[M_{a,j}^p]<∞ whenever a<2 and 1≤p<4/a². (2)

View M_{a,j} as the integral over the unit square of the exponential of the piecewise constant Gaussian field X_j(z)=h_{S_j(z)} with its pointwise variance subtracted. Uniformly in j,z,w,

    Cov(X_j(z),X_j(w))≤log(1/(|z−w|∨2⁻ʲ))+C₀.  (3)

Construct a star-scale Gaussian field Y_t on a neighborhood of the closed square, with independent scale increments and covariance

    Cov(Y_t(z),Y_t(w))=∫₀ᵗ exp(−e^{2s}|z−w|²/2)ds.

The integrand is positive definite in z−w. At t=jL its covariance equals log(1/(|z−w|∨2⁻ʲ))+O(1) uniformly on the square. Add one independent centered scalar Gaussian of sufficiently large fixed variance to Y_{jL}; its covariance now dominates (3) pointwise. Kahane's convexity inequality (Rhodes–Vargas, review, Theorem 2.1), applied to the convex function u↦u^p, compares E[M_{a,j}^p] to the pth moment of the star-scale total mass times a fixed lognormal factor. The star-scale masses form a positive martingale. Rhodes–Vargas Theorem 2.5 gives its uniform integrability for a<2, and Theorem 2.11 gives a finite pth moment of its limit for p<4/a². Uniform integrability identifies each mass as the conditional expectation of the limit, so conditional Jensen bounds all its pth moments by that of the limit. This proves (2).

Now fix any p<x. Since γ<2, we can choose a>γ sufficiently close to γ that a<2 and p<4/a². The quadratic β(a)=2−aQ+a²/2 has roots γ and 4/γ, so β(a)<0. Minkowski's inequality, (1), and (2) imply

    ||I||_p≤Σ_j||I_j||_p≤CΣ_j exp(β(a)jL)<∞.

Thus E[N^p]<∞. Markov's inequality gives P(N≥n)≤E[N^p]n⁻ᵖ. Letting p increase to x proves the logarithmic upper bound. ∎

**Sources and scope.** APPS Definition 6.2 specifies the active condition: https://arxiv.org/html/2005.11845. Rhodes–Vargas, *Gaussian multiplicative chaos and applications: a review*, Theorems 2.1, 2.5, and 2.11: https://arxiv.org/html/1305.6221. The argument only needs the local logarithmic Green kernel and gives no lower tail. It also proves the earlier second-moment theorem as a special case when x>2.
