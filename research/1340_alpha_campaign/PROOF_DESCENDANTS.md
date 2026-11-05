# High surplus creates exponentially many reachable descendants

Status: proved; independently reviewed by alpha_entropy, with the ball-domain clarification below incorporated. This is a conditional estimate for a deterministic square and its own average, not for a square chosen by an exploration. The marking argument in section 5 explains how to use it without assuming that stronger conditional law.

## 1. Statement

Use the curvature-normalized APPS sphere GFF, fixed threshold epsilon>0, and Q>2. Write L=log2 and gamma=Q−sqrt(Q²−4). Fix a>Q. For a deterministic dyadic square S of depth j, condition on h_S=h with h∈[ajL,ajL+1]. Put

    y0=(a−Q)jL−Q log epsilon.

For every eta>0, uniformly over S and such h, the conditional probability tends to one as j→infinity that at least exp((gamma−eta)y0) descendant squares of S are active along their entire paths starting at S. In particular, if S is reachable in the original subdivision, all these squares are split there. No reachability assumption is included in the conditioning.

More generally, the same argument works for Q≥2, a>Q and any fixed b∈(0,2) with b<Q. In that form replace gamma in the conclusion by

    R_b=(2−b²/2)/(Q−b).

This extension includes Q=2, when R_b=(2+b)/2. The Gaussian field at this endpoint is well defined; the proof uses only auxiliary subcritical b-chaos, not critical GMC. This generalized statement passed the separate independent review in `1343_alpha_meander_bridge/REFEREE_CRITICAL_LOWER.md`. Its requested clarification about the ceiling-restricted second-moment denominator is incorporated in section 5.

## 2. Uniform conditional field estimates

The sphere Green kernel on a neighborhood of the unit square is −log|z−w|+H(z,w), with H bounded smooth. Average it over dyadic squares. For descendant U⊂S at relative depth s, and descendants U,V at equal relative depth ell, logarithmic rescaling gives, uniformly in all squares,

    Var(h_S)=jL+O(1),
    Cov(h_U,h_S)=jL+O(1),
    Var(h_U)=(j+s)L+O(1).

For nested descendants U_s⊃U_ell, Cov(h_(U_s),h_(U_ell))=(j+s)L+O(1). For arbitrary U,V at relative depth ell, after scaling S to the unit square their covariance is bounded above by

    jL+log(1/(|z−w|∨2^(−ell)))+C

for z∈scaled U,w∈scaled V. For adjacent or touching squares, these bounds follow by scaling the integral of −log distance; its integral against a unit-square uniform measure is uniformly bounded, even when the other square shrinks to a point in its closure. For separated squares, distance is comparable to the center separation, except within a bounded number of neighboring cells where the preceding argument applies. Smooth H contributes a uniform bound.

Let X_U=h_U−E[h_U|h_S=h]. Gaussian regression and h=ajL+O(1) give

    E[h_U|h_S=h]=h+O(1),
    Var(X_U)=sL+O(1),
    Cov(X_(U_s),X_(U_ell))=sL+O(1),              (1)

and the conditional equal-depth covariance, in the scaled coordinates, is bounded above by log(1/(|z−w|∨2^(−ell)))+C. All constants are uniform in j,S,h. Indeed the rank-one regression subtraction is jL+O(1). For sufficiently large j its denominator Var(h_S) is positive and comparable to jL. The conditional covariance does not depend on h.

## 3. Negative moments of the scale mass

Fix b∈(0,2), eventually b=gamma. Define

    w_U=4^(−ell)exp(bX_U−b²Var(X_U)/2),
    M_ell=Σ_(U at relative depth ell) w_U.

For each fixed r>0 there is C_r<infinity such that

    E[M_ell^(−r)|h_S=h]≤C_r                    (2)

uniformly in j,S,h,ell. Here is a justification that avoids assuming the square-average approximation is a martingale. Regard X_U as a piecewise constant field on the scaled unit square. Compare it with the independent-scale field Y_t having covariance

    K_t(z,w)=integral_0^t exp(−e^(2s)|z−w|²/2) ds,

at t=ell L, plus an independent scalar Gaussian of fixed sufficiently large variance. K_t equals log(1/(|z−w|∨e^(−t)))+O(1) uniformly on the square. Thus it dominates the covariance in section 2 after this scalar addition. Kahane's inequality for the convex decreasing function u↦(u+v)^(−r), followed by v↓0, bounds the left side of (2) by the corresponding negative moment for Y_t with its independent lognormal scalar factor.

Define Y_t on an open neighborhood of the closed square. Its normalized total masses form a nonnegative martingale in t. Since b<2, the subcritical nondegeneracy and uniform-integrability theorem makes this martingale the conditional expectation of its limiting mass. The limiting square mass is at least the mass of a fixed inscribed compact ball, whose negative moments of every order are finite. Conditional Jensen therefore bounds each truncated square-mass negative moment by the limiting one. Adding the scalar factor multiplies the bound by a fixed finite lognormal moment. These are Rhodes–Vargas, *Gaussian multiplicative chaos and applications: a review*, Theorems 2.1, 2.5 and 2.12, https://arxiv.org/html/1305.6221 .

For completeness the integral version of Kahane follows from its finite Gaussian-vector form: partition each of the finitely many square cells into a refining grid, use point values of the smooth comparison field and the constant value of X on each cell, and pass to integrals. At fixed t the comparison field is smooth on a compact neighborhood and its supremum has Gaussian tails, which justifies convergence of the comparison negative moments. The left mass is unchanged by this refinement. Regularization by v>0 also makes the passage bounded before v↓0.

## 4. Discarding bad paths loses exponentially little mass

Choose fixed b∈(0,2) with b<Q (b=gamma for the initial subcritical conclusion), and fixed small delta,zeta>0. Take

    ell=floor((1−delta)y0/((Q−b)L)).

A terminal descendant is good if every intermediate ancestor at relative depths 1,...,ell is active and its terminal centered value satisfies X_U≤(b+zeta)ell L. For each terminal U, the weighted Gaussian law with density exp(bX_U−b²Var(X_U)/2) shifts the mean of X_(U_s) to b Cov(X_(U_s),X_U)=bsL+O(1), leaving its variance unchanged. By (1), inactivity at depth s requires

    X_(U_s)≤QsL−y0+O(1).

The difference between tilted mean and this threshold is at least

    y0−(Q−b)sL−O(1)≥delta y0−O(1).

Since all the variances are at most ell L+C=O(y0), a Gaussian tail bound and a union bound over s give a bound exp(−c y0) on failure of a path, absorbing the O(y0) number of depths. The tilted terminal mean is b ell L+O(1); failure of the terminal upper bound has the same exponential estimate. Consequently, if B_ell is the total weight of bad terminal descendants,

    E[B_ell|h_S=h]≤exp(−c y0)                  (3)

for some c=c(Q,a,delta,zeta)>0 and all sufficiently large j, uniformly over S and h.

Choose theta>0 with theta<c/2. From (2) and Markov,

    P(M_ell<e^(−theta y0)|h_S=h)≤C_r e^(−r theta y0),
    P(B_ell>(1/2)e^(−theta y0)|h_S=h)
       ≤2 exp(−(c−theta)y0).

Thus with conditional probability tending to one, the good mass is at least (1/2)e^(−theta y0). Each good terminal weight is at most

    C exp(−(2−b²/2−b zeta)ell L).

The number of good terminal descendants is therefore at least

    (1/(2C)) exp((2−b²/2−b zeta)ell L−theta y0). (4)

The exponent in (4), divided by y0, tends to R_b as delta,zeta,theta↓0. When Q>2 and b=gamma, the relation 2−gamma²/2=gamma(Q−gamma) identifies R_b with gamma. For any given eta>0, first choose delta,zeta sufficiently small, then choose theta>0 smaller than both c/2 and the remaining allowed error. This proves section 1, including its generalized form. The probability convergence is uniform, even though its rate depends on eta.

## 5. Correct way to attach this to a high-spine second moment

Suppose a high-spine proof counts deterministic depth-j squares satisfying h_S∈[ajL,ajL+1] and ancestral reachability. Suppose it supplies a uniform eta0>0 for the probability of ancestral reachability conditional on the square's high band. Mark a candidate as successful when its descendants also satisfy section 1. Uniform conditional success tending to one gives

    P(reachable and successful | high band)≥eta0−o(1).

Hence the first moment of the successful count retains the high-band first-moment exponent. A second-moment bound must be proved for the precise dominating count used; existence of a high-spine event alone does not suffice. The current general proof in `1341_alpha_tail/PROOF_spine_hit.md` additionally restricts each candidate by an upper ancestral ceiling. Conditional on the endpoint high band this ceiling fails with probability o(1), so the reachable, ceiling-restricted, descendant-successful count still has expectation at least (eta0−o(1)) E[Z_j], where Z_j is the all-high-band count. Its second moment is bounded by that of the ceiling-restricted count T_j, for which the proof gives

    E[T_j²]≤poly(j) exp(a delta j log2) E[Z_j].

The unrestricted count Z_j has the stronger bound without the exponential penalty only in the earlier a²≥8 regime; it must not be used that way for general a>2. Paley–Zygmund with T_j and then delta↓0 yields the required exponential lower bound. This bypasses the unjustified conditional law given a first-found square and its full exploration history.

For a=rho_plus=Q+sqrt(Q²−4), section 1 supplies exp((gamma−eta)y0) descendants. Combined with a high-spine hit probability exp(−a y0−o(y0)) obtained by this marked-count method, it implies the target lower tail exponent a/gamma=4/gamma². This last conclusion awaits the independent high-spine proof and its uniform ancestor-survival check.
