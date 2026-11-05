# The continuum exponent curve, including the critical endpoint

Status: proved, with separate independent agent reviews of the ingredients and assembly. This is an internal research proof, not a claim of external peer review or established publication novelty.

## Model and tail theorem

Fix a smooth metric g0 on the Riemann sphere that is flat on a neighborhood of the closed unit square. Let h be its Gaussian free field with Dirichlet inner product (2pi)^(-1) integral grad f dot grad g, normalized by (h,K0)=0, where K0 is Gaussian curvature. For a dyadic square S, let |S| be its side length and h_S its area average. Fix epsilon>0 and Q≥2. Starting from the unit square, split a square into four children when

    exp(h_S/Q)|S|>epsilon.

Let N_Q be the terminal-leaf count. Put gamma=Q−sqrt(Q²−4), taking gamma=2 at Q=2. Then

    lim_(n→infinity) log P(N_Q≥n)/log n = −4/gamma².       (T)

The limit is at fixed Q, epsilon and background metric. It asserts logarithmic asymptotics, not a nonzero limiting constant. The event is well defined for all integers n; exact leaf counts lie on 1 modulo 3. At Q=2, E[N_Q]<infinity as well, so the theorem cannot be strengthened to a pure constant-times-1/n tail without additional logarithmic corrections.

The subdivision and field are exactly APPS Definitions 6.2 and 6.7, https://arxiv.org/html/2005.11845 . This is the correlated GFF model, not a branching cascade.

## Determinant theorem

On N_Q=m, let h^epsilon be the finite conditional projection from APPS Definition 6.4 and set

    V=Vol_(exp(2h^epsilon/Q)g0)(sphere),
    C=(1/4)log m−(3/4)log V,
    g=exp(2(h^epsilon/Q+C))g0,
    A0=det'_zeta Delta_g0 / Vol_g0(sphere).

Define det'_PA Delta_g using the Polyakov–Alvarez formula (APPS equation 6.7). This avoids assuming spectral-zeta continuation for the nonsmooth random metric. In orthonormal projection coordinates x_S=alpha_S/Q,

    det'_PA Delta_g=A0 m^(1/3) exp(−sum_S x_S²/6).        (D)

At base Q0²=25/6 (central charge zero), for every real k≥−1/2 set Q_k²=Q0²+k/3 and define

    C_k(n)=sum_(m≥n) (Q_k/Q0)^m A0^(−k)m^(−k/3)
                     E_(Q0)[(det'_PA Delta_g)^k;N_(Q0)=m].

The finite-dimensional Gaussian tilt gives the exact identity

    C_k(n)=P_(Q_k)(N_(Q_k)≥n).                            (G)

Consequently

    C_k(n)=n^(−x_k+o(1)),
    x_k+1/x_k=(13+2k)/6,   x_k≥1,
    x_k=(13+2k+sqrt((13+2k)²−144))/12.                   (K)

At k=2, x_2=(17+sqrt145)/12. Thus 1+x_2=(29+sqrt145)/12 is the number predicted for the meander exponent, now obtained from a cumulative partition-function theorem in this specified continuum model. Equation (K) does not by itself identify any exponent of discrete meander counts.

The compensation in C_k is part of the theorem. For the energy-normalized determinant D=exp(−sum x_S²/6), it reduces to the exponential compensation in folder 159. For the actual PA-defined determinant, omitting m^(−k/3) changes the cumulative exponent to x_k−k/3. The distinction is invisible when conditioning on one fixed m but essential across sizes.

## Proof structure and independent checks

| Step | Argument | Proof | Independent review |
|---|---|---|---|
| Subcritical upper | Dominate all active squares at each scale by decaying multiples of auxiliary subcritical chaos; positive moments and Minkowski | 1341/PROOF_sharp_upper.md | 1343/REFEREE_SHARP_UPPER.md |
| Reachable high square | High endpoint band plus upper ancestral ceiling; Euclidean pair shells and a three-average Gaussian Chernoff bound | 1341/PROOF_spine_hit.md | 1343/REFEREE_SPINE_GENERAL.md |
| Descendant growth | Negative chaos moments and tilted bad-path estimates give conditional success tending to one at a deterministic high square | 1340/PROOF_DESCENDANTS.md | 1342/REFEREE_APPS_DESCENDANTS.md; 1343/REFEREE_CRITICAL_LOWER.md |
| Subcritical assembly | Mark deterministic candidates; retain the first moment and use the ceiling-restricted second moment | 1341/PROOF_sharp_tail_all_subcritical.md | 1340/REFEREE_SHARP_TAIL_ASSEMBLY.md |
| Critical upper | Endpoint Gaussian tilt, square-to-circle comparison and Brownian reflection yield summable expected reachable-depth counts | 1342/PROOF_APPS_CRITICAL_UPPER.md | 1340/REFEREE_CRITICAL_UPPER.md |
| Critical lower | High-band slope a>2 and auxiliary chaos tilt b<2, then a decreases to 2 and b increases to 2 | 1341/PROOF_critical_lower.md | 1343/REFEREE_CRITICAL_LOWER.md |
| Endpoint assembly | Match logarithmic upper and lower bounds | 1341/PROOF_critical_tail.md | 1340/REFEREE_CRITICAL_ASSEMBLY.md |
| Exact tilt | Curvature normalization, full rank m Gaussian density, and explicit polynomial determinant factor | 1340/PROOF_NORMALIZATION.md | 1343/REFEREE_NORMALIZATION.md |

Folder abbreviations in the table mean 1340_alpha_campaign, 1341_alpha_tail, 1342_alpha_tree_entropy, and 1343_alpha_meander_bridge.

The upper-bound chaos inputs are Kahane comparison, subcritical uniform integrability, and positive moments; the descendant lemma also uses negative moments. These are Rhodes–Vargas, *Gaussian multiplicative chaos and applications: a review*, Theorems 2.1, 2.5, 2.11 and 2.12, https://arxiv.org/html/1305.6221 . The critical proof does not invoke a critical GMC measure. All probability estimates use the actual GFF covariance; no independent-branch approximation is used.

## Claims not supplied

No local law P(N=m)=m^(−1−x+o(1)) on the attainable lattice is proved. No interchange of the k derivative with an asymptotic error is proved. Thus the discrete −4/9 log n correction does not follow by differentiating (K). No convergence theorem from determinant-weighted maps or meanders to this model is proved. The BGS unit-area decorated scaling conjecture alone does not determine counting normalization. Those are separate open steps.
