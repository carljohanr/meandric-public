# Exact balanced-area exponent in a conditional two-disk model

**Status.** Independently reviewed in `1342_alpha_tree_entropy/REFEREE_TWO_DISK_BALANCE.md`. Complete calculation for two independent standard quantum-disk areas at a prescribed common boundary length, conditioned on their total area. It is a mechanism for the proposed continuum exponent, **not** a proof of the Borga–Gwynne–Sun meandric-permuton small-ball law. No identification of the off-block event with this disk model is asserted.

Put `a=4/γ²` and fix `γ∈(0,2)`; this exponent `a` is unrelated to Ang–Gwynne's variance constant denoted by the same letter. Ang–Gwynne, *Liouville quantum gravity surfaces with boundary as matings of trees*, arXiv:1903.09120, Theorem 1.2, gives the area `T` of the unit-boundary-length quantum disk the inverse-gamma density

    f₁(t)=K t^(−1−a) exp(−β/t),     t>0,                 (1)

for positive constants `β` and `K=β^a/Γ(a)` (their displayed constants give this normalization). The exact scaling of LQG area and boundary length under a field constant shift gives the density for boundary length `ℓ`:

    f_ℓ(t)=K ℓ^(2a)t^(−1−a)exp(−βℓ²/t).              (2)

Let `A₁,A₂` be independent with density `f_ℓ`, and disintegrate their joint law at `A₁+A₂=1`. This is a well-defined conditional density because the convolution `g_ℓ(1)=(f_ℓ*f_ℓ)(1)` is positive and finite. For `U=A₁ | (A₁+A₂=1)`,

    p_ℓ(u)=f_ℓ(u)f_ℓ(1−u)/g_ℓ(1),  0<u<1.           (3)

The denominator has the sharp asymptotic

    g_ℓ(1) ~ 2K ℓ^(2a).                              (4)

To see this, fix `η∈(0,1/2)`. On `u∈[η,1−η]`, both densities are `O(ℓ^(2a))`, so that contribution is `O(ℓ^(4a))`. For `u∈(0,η)`, `f_ℓ(1−u)/(Kℓ^(2a))` tends uniformly as `ℓ→0` to `(1−u)^(−1−a)` and is bounded by a constant depending on `η`. Meanwhile `A₁=ℓ²T` tends to zero in probability, so the left endpoint contribution divided by `Kℓ^(2a)` tends to one. The right endpoint contributes the same. Thus (4) follows.

For fixed `v∈(0,1)` and `δ→0` jointly with `ℓ→0`, provided `δ=o(1)` and the window `[v−δ,v+δ]` stays inside `(0,1)`, (2)–(4) and uniform continuity on a compact interval yield

    P(|U−v|≤δ)
      ~ K [v(1−v)]^(−1−a) δ ℓ^(2a).                  (5)

The factor `2δ` from the window cancels the factor `2` in (4). At `v=1/2` and the proposed neck-resolution relation `ℓ²=δ`, equation (5) becomes

    P(|U−1/2|≤δ | A₁+A₂=1)
      ~ K 4^(1+a) δ^(1+a),       a=4/γ².             (6)

More generally `ℓ²` between fixed positive multiples of `δ` gives matching constant-factor bounds and the same logarithmic exponent `1+4/γ²`.

This calculation also distinguishes two costs. Given the total-area-one condition and a small boundary length, one disk is typically microscopic and the other carries nearly all area; the relative cost of *two* macroscopic disks is `ℓ^(2a)`. Restricting their split to a window of width `δ` around one half then costs another factor `δ`. The heavy area tail from Theorem 1.2 supplies the first cost, while a continuous area-split density on compact subintervals supplies the second. A naive multiplication of two unconditioned disk tails would instead miss the normalization of the total-area-one disintegration.

**Missing geometric bridge.** Borga–Gwynne–Sun, arXiv:2207.02319, Definition 1.7, defines the candidate meandric permuton using a unit-area `γ`-quantum sphere and two independent space-filling `SLE₈` curves, with `γ²=(17−√145)/3`. The relevant event `E_δ` is a near block in their two area visit orders. To apply (6), one would need all of the following, with logarithmic-relative errors: (i) an area-resolution-`δ` geometric neck whose quantum boundary length is comparable to `√δ`; (ii) a sphere neck disintegration comparable to two independent *standard* length-`ℓ` disks conditioned on total area one, including its boundary-length selection measure; and (iii) a conditional SLE₈ estimate that both independent curves traverse the two blobs in block order with probability `δ^{o(1)}` and that off-block mass is at most `δ`. None is supplied by Theorem 1.2.

The available welding and mating theorems have different curve parameters. Ang–Holden–Sun, arXiv:2009.08389, Theorems 2.2–2.4, weld disks along their quantum boundary via `SLE_{γ²}`-type interfaces; their Theorem 7.5 gives a sphere Brownian boundary-length encoding for a matched space-filling curve. Ang–Gwynne Theorems 1.1–1.3 likewise describe the matched space-filling curve `SLE_{16/γ²}`. At the meandric value, `16/γ²≈9.68`, not 8. These theorems therefore do not identify the law of the two independent `SLE₈` visit orders or prove (i)–(iii). In particular the value of `γ` alone does not determine the small-ball exponent of `E_δ`.

There is not even a deterministic inclusion from a block-order event to a thin geometric neck. On an ordinary smooth sphere, divide the sphere into two hemispheres with a macroscopic equatorial boundary. Choose two continuous area-parametrized space-filling curves that each fill the first hemisphere during time `[0,1/2]` and the second during `[1/2,1]`. Their induced permuton has zero off-block mass for the half-interval aligned with these times, although the separating equator has fixed positive length. This construction is **not** an SLE₈ sample; it shows that a thin-neck conclusion requires a probabilistic theorem about the SLE₈-decorated law, rather than the definition of `E_δ` alone.
