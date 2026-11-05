# Exact reference-circle change for the BGS sphere law

Status: independently reviewed in `1342_alpha_tree_entropy/REFEREE_SPHERE_PINNED_CLUMPS.md`, against BGS Definition 2.2. This resolves a normalization issue in the proposed two-clump lower bound, not the whole small-ball theorem.

Write a=4/gamma^2 and beta=a-2. Let h_1 be the whole-plane GFF with average zero on the unit circle, as in BGS (2.2)–(2.3), and take R>3. Set h_R=h_1-(h_1)_R and c=(h_R)_1=-(h_1)_R. Then h_1=h_R-c exactly. The covariance of h_R is

    G_R(z,w)=log((|z| vee R)(|w| vee R)/(R|z-w|)).

The scalar c is centered Gaussian with variance log R. By circle averaging the logarithmic kernel,

    Cov(h_R(z),c)=log(|z| vee R)-log(|z| vee 1).

These statements follow either by subtracting the R-circle average twice from the unit-circle covariance or by the elementary logarithmic potential of uniform circle measure.

Let d be the deterministic part of the BGS Liouville field. Expanding its definition gives

    d(z)=-gamma log|z|-gamma log|z-1|-beta gamma log(|z| vee 1).

Indeed 3gamma-2Q=-beta gamma. If M_R is the total area of h_R+d, then the total area of h_1+d is exp(-gamma c)M_R. Normalized measures are identical, since this common multiplier cancels. Thus every nonnegative functional F of normalized area satisfies the exact sphere expectation identity

    E_sphere F = Z^(-1) E_R[exp(-beta gamma c) M_R^beta F],

where Z is the positive finite BGS normalizer. Here LQG area uses the field's circle-average renormalization, so adding a constant multiplies the area by its exponential exactly.

Apply the elementary Gaussian exponential-tilt identity to the joint Gaussian field h_R and c. Weight exp(-beta gamma c) changes the mean of h_R by

    m(z)=-beta gamma[log(|z| vee R)-log(|z| vee 1)]

and contributes the fixed constant exp(beta^2 gamma^2 log R/2). Hence the same normalized sphere law is obtained by taking a centered h_R, adding

    d_R(z)=d(z)+m(z)
          =-gamma log|z|-gamma log|z-1|-beta gamma log(|z| vee R),

and weighting its total area by the beta-th power, with the corresponding fixed normalizing constant. This identity is for arbitrary nonnegative measurable F, by Gaussian tilting and monotone approximation; it is not a comparison that loses an exponent.

Inside |z|<R, the final term is the constant -beta gamma log R. There are just two gamma insertions and a smooth constant background near 0 and 1. Outside R the decay term remains; it must not be discarded when defining the global area and sphere weight.

The reference circle now lies outside any two sufficiently small disjoint disks centered at 0 and 1. Their zero-boundary Dirichlet components have zero trace on the R-circle and are independent of the exterior in the usual GFF domain Markov decomposition. This is the independence needed by the two-clump argument. Using the original unit-circle representative directly would not justify that assertion for the disk centered at 1, since it intersects the normalization circle.

## Conditional annulus control without boundary regularity assumptions

For a compact annulus away from both insertions, choose s with a-1<s<a. The unconditioned s-th area moment is finite by subcritical GMC moments. Consequently the conditional s-th moment given the exterior is finite almost surely, by conditional expectation. This remains true for the residual field after extracting an independent bounded smooth Gaussian coordinate: on the event |Z|<=1, its GMC mass and the original mass differ by bounded factors; independence of Z and the residual converts finiteness of the original moment into finiteness of the residual moment.

Thus for almost every exterior realization, conditional Markov's inequality bounds the annulus mass above epsilon t by O_ext(t^-s)=o(t^(-(a-1))). There is no need to demand bounded harmonic boundary data up to the Dirichlet disk boundary. On each strictly smaller core the harmonic term is smooth, so a pointwise positive pinned-band asymptotic can be combined with this annulus estimate. Fatou's lemma over a positive-probability exterior event yields a joint lower bound without a uniform family version of the pinned-tail theorem.

This last step assumes the local pinned-band theorem has been separately verified. It does not replace its proof or establish the final Gaussian balance event.

Primary source: Borga–Gwynne–Sun, arXiv:2207.02319v3, Definition 2.2, equations (2.2), (2.3), (2.9), and (2.10): https://arxiv.org/html/2207.02319v3 .
