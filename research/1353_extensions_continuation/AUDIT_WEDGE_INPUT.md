# Original-source check and focused audit of the wedge input

Status: original source verified and the required local estimate rederived. Independent review requested. This closes the specific source-provenance gap recorded by 639/650; it is not a second full review of the entire band-exponent proof.

## Source and normalization

Bañuelos and Smits, Brownian motion in cones, Probability Theory and Related Fields 108 (1997), 299–319. Author manuscript: https://www.math.purdue.edu/~banuelos/Papers/cones.pdf . Lemma 1, equation (2.2), is on manuscript page 6 (PDF page 6). Lemma 2 specifies the generator as one half of the Laplacian. For a planar wedge of angle theta, the angular eigenfunctions are sqrt(2/theta) sin(j*pi*phi/theta), and the angular eigenvalues are (j*pi/theta)^2. Substitution in that equation gives exactly the factor 2/(theta*t), Bessel order j*pi/theta and Lebesgue-area kernel used in 639 REPAIRS equation (2.1). No polar Jacobian belongs in the kernel itself; integration uses r' dr' dphi'.

The source was read directly, not inferred from numerical agreement or a secondary quotation. This verifies the literature input left unchecked in 639 FINAL and inherited in 650 FINAL. We do not need a separate citation to Carslaw–Jaeger or its section numbering.

## The local estimate actually needed

Use the notation and proof strategy of 639 REPAIRS, Proposition B2. Let V be the open wedge, nu=pi/theta, let Z,D be compact subsets of V with D having positive area, and fix 0<t0<t1<infinity. There is c>0 such that

P_(lambda*z)(tau_V>t, B_t in D) >= c lambda^nu

for z in Z, t in [t0,t1], and 0<lambda<=1. Here tau_V is the first exit time, and B is standard planar Brownian motion. In this statement compact subsets of V stay a positive distance from its sides and apex.

Proof of the specialization. Write x=(r,phi), y=(s,psi) in polar coordinates and w=rs/t. All radii of z and y stay in fixed positive compact intervals and all angles stay away from the wedge sides. The first angular term has sin(nu*phi)sin(nu*psi)>=c0>0. The elementary Bessel series gives I_nu(w)>=(w/2)^nu/Gamma(nu+1). The absolute value of all higher angular terms is at most
 e^w sum_(j>=2)(w/2)^(j*nu)/Gamma(j*nu+1).
For 0<w<=w0, this sum is O(w^(2nu)), uniformly: factor w^(2nu) and bound the remaining convergent Gamma-denominator series at fixed w0. Since w is comparable to r on the compact target/time ranges, the first term dominates for sufficiently small r. The common Gaussian and 1/t factors are bounded below there. Thus p_t^V(x,y)>=c r^nu uniformly for y in D and small r. Integrate over D. For the remaining r range, positivity and continuity of the killed heat kernel on an interior compact set give a positive minimum. This proves the displayed bound; no fitted exponent or simulation is involved.

This proof explicitly controls the signs of higher angular terms rather than assuming every term of the sine series is positive. The estimate is uniform in both the small radius and the time in its compact interval.

## The stage constant required by 643/650

The relevant covariance matrix Sigma is positive definite. An invertible linear map turns its Brownian motion in the quadrant into standard Brownian motion in V. Compact interior boxes remain compact interior sets, and the constant Jacobian changes only multiplicative constants.

Let delta in (0,1/8), time u in [3/4,1], initial point zeta in [1,4]^2, and margin v with 0<=v<=sqrt(delta)/8. Subtract v(1,1). The starting points lie in a fixed interior compact set. The translated target contains sqrt(delta)[5/4,3]^2. The required path event is containment in the quadrant and arrival in that target.

For a fixed interior intermediate compact set D', symmetry and the semigroup identity imply
 p_u(zeta',y) >= inf_(w in D') p_(u/2)(w,zeta') * P_y(tau>u/2, W_(u/2) in D').
The first factor has a fixed positive lower bound over the start and time ranges. The local estimate just proved bounds the second below by c delta^(nu/2) for y in sqrt(delta)[5/4,3]^2. This target has area proportional to delta. Integration therefore yields c4 delta^(1+nu/2), with c4 independent of delta, u and margin v. Since the older notation is mu_wedge=nu/2, this is exactly the sharp stage input BM4 and St1. Boundary contact conventions can be handled by taking strict containment; Brownian endpoint boundary sets have zero measure.

The factor c4 is chosen before delta in the multiscale lower proof. This independence is the key reason the per-stage constant contributes only log(1/c4)/log(1/delta), which can be made small. The source check supplies that required uniformity, not merely a fixed-delta lower bound.

## Other inherited analytic points

Whitt, Proofs of the martingale FCLT, Probability Surveys 4 (2007), Theorem 2.1(ii), https://arxiv.org/pdf/0712.1929 , verifies the martingale theorem used by 639. The existing 639/650 argument time-changes by the first move reaching a point-time; bounded increments and predictable covariation convergence are its substantive hypotheses. Their proofs and the corresponding algebraic identities remain inherited reviewed inputs. Reading the theorem alone does not reverify those identities or the uniform regular-state argument.

The lower exponent further uses: the cylinder version of reversal, a factor two for reversed stopping times, overshoot in the conditioned state, the polynomial fallback, and Jensen only after finite logarithms are established. These repairs are recorded explicitly in 650 Sections 1–5. The new corollary must cite 643 together with 650, not the original uncorrected argument by itself.

No claim of a Brownian conditioned excursion limit is used or justified by these estimates. E1 only needs the two-sided logarithmic coefficient exponent. The sharp constant, full asymptotic equivalent, and band permuton convergence are separate questions.
