# One-sided small-ball lower bound for the BGS meandric permuton

**Status.** Proved and independently reviewed. The pinned-GMC local tail in `1342_alpha_tree_entropy/PROOF_PINNED_GMC_LOCAL_TAIL.md` was reviewed in `1343_alpha_meander_bridge/REFEREE_PINNED_GMC_LOCAL_TAIL.md`; the sphere gauge and clump assembly were reviewed in `1342_alpha_tree_entropy/REFEREE_SPHERE_PINNED_CLUMPS.md`; the fixed-clump SLE implication was reviewed in `1342_alpha_tree_entropy/REFEREE_FIXED_CLUMP_SLE.md`. No upper bound or exact continuum exponent is claimed here.

Let `γ∈(0,√2)` be the Borga–Gwynne–Sun (BGS) sphere parameter, `a=4/γ²>2`, and `E_δ` the trimmed off-block event of `PROOF_OFFBLOCK_SET_IDENTITY.md`. There are constants `c,δ₀>0` depending on `γ` such that

`P(E_δ)≥c δ^(1+a)` for all `0<δ<δ₀`.       (1)

Consequently,

`limsup_(δ↓0) log P(E_δ)/log δ ≤ 1+a`.       (2)

At the meandric-permuton parameter `γ²=(17−√145)/3`, the right-hand side exponent is `(29+√145)/12`. This is a continuum probability bound, not a finite-meander statement.

## 1. Exact outer-circle gauge for the triply marked sphere

BGS Definition 2.2 gives the unnormalized Liouville field

`h_L=h^C+γG_C(0,·)+γG_C(1,·)−(2Q−γ)log|·|_+`, `Q=2/γ+γ/2`,

where `h^C` has zero average on the unit circle. The unit-area law normalizes its GMC by total mass `M` and biases the underlying field law by `M^κ`, `κ=a−2>0`. Choose `R>3`. Let `h_R` be the same whole-plane GFF modulo constants anchored to have zero average on the circle of radius `R`, and `c=(h_R)_1`. Then `h^C=h_R−c`. For any nonnegative functional `F` of the normalized area measure,

`E_sphere[F]=Z^(-1)E_R[e^(−κγc) M(h_R+d)^κ F(normalize M(h_R+d))]`,

with `d=γG_C(0,·)+γG_C(1,·)−(2Q−γ)log|·|_+`. Gaussian tilting by `e^(−κγc)` changes `h_R` to `h_R+m`, where `m(z)=−κγ Cov(h_R(z),c)`. The additive Gaussian normalizer cancels from `Z`. The circle potential gives `Cov(h_R(z),c)=log R−log|z|_+` for `|z|<R` and zero for `|z|≥R`. Since `κγ=4/γ−2γ`, direct expansion yields the exact deterministic background

`d(z)+m(z)=−γlog|z|−γlog|z−1|−κγlog(|z|∨R)`.       (G)

Thus under the sphere law we may use the `R`-anchored GFF plus (G), bias its unnormalized total mass to power `κ`, and then normalize area. Equation (G) has clean `γ` insertions at `0` and `1`; its last term is constant in a neighborhood of both. It must be retained outside `R` when defining the global total mass.

## 2. Two independent pinned local masses

Fix disjoint balls `B_i=B_r(i)` with closures in `B_R(0)`, and core balls `K_i=B_(r/2)(i)`. The `R`-anchored GFF has the domain Markov property in the balls: conditional on the exterior field, its restrictions are independent zero-boundary GFFs plus harmonic functions. Write `A_i=M(B_i)`, `O=M(C\(B_0∪B_1))`. The exterior mass `O` is finite almost surely. Choose an exterior event of positive probability on which `O≤B` and the harmonic functions are bounded on the cores. Constants may depend on that event and on fixed `r,R`.

Extract from the first zero-boundary GFF a scalar `Z_0` along a smooth nonnegative Cameron–Martin function `f`, supported inside `B_0` and equal to one on `K_0`. Its Gaussian density is bounded below on every fixed compact interval. The residual first field and the entire second-ball field are independent conditional on the exterior. Include the deterministic GMC Wick factor for `Z_0` in the first residual core mass `C_0`, so that at scalar value `z` the actual first-core mass is exactly `e^(γz)C_0`.

The pinned-GMC local tail theorem in `1342_alpha_tree_entropy/PROOF_PINNED_GMC_LOCAL_TAIL.md` applies pointwise for almost every exterior field to both core masses (the first with a smooth rank-one covariance subtraction). With `q=a−1`, it gives positive finite conditional constants `c_i` such that

`P(t≤C_i≤2t | exterior)~c_i t^(−q)`.       (3)

The annuli `B_i\K_i` stay a fixed positive distance from the relevant insertion. Their masses have finite conditional `s` moments for some `q<s<a`, because unconditional moments are finite and Fubini applies. The same holds for the residual first annulus, since on the independent event `|Z_0|≤1` the full and residual annulus masses differ by a bounded factor. Thus each annulus exceeds `t/10` with conditional probability `o(t^(−q))`. Intersecting this with (3) leaves the same leading order. No uniform local-tail theorem over exterior fields is required: restrict further to a positive-probability exterior subset where `c_0,c_1` are bounded below and the conditional annulus moments bounded above.

Conditional independence now shows that the event `C_i∈[t,2t]`, both annuli≤`t/10`, and `O≤B` has probability at least `c t^(−2q)` under the unbiased outer-gauge field law. More explicitly, let `g_t(exterior)` be `t^(2q)` times the product of the two conditional good-band probabilities. For almost every exterior on the selected set, `liminf g_t≥c_0 c_1>0`. Fatou's lemma gives `liminf_(t→∞) E[1_selected g_t]≥E[1_selected c_0 c_1]>0`; hence one fixed positive constant works for every sufficiently large `t`. This avoids any implicit uniformity over random harmonic boundary fields.

## 3. Gaussian balance coordinate and sphere bias

Condition on this event and vary `Z_0=z`. The first-ball mass `A_0(z)` is continuous and increasing. Its core contribution is `e^(γz)C_0`; the annulus contribution stays `O_K(t)` for `z∈[−K,K]`. Choose a fixed sufficiently large `K` so that, uniformly over the bulk bands and annulus bounds,

`A_0(−K)<A_1+O<A_0(K)`.

Also `A_0'(z)≤C_K t` on that interval, since `0≤f≤1`. Therefore the set of `z` satisfying `|A_0(z)−A_1−O|≤1` contains an interval of length at least `c/t`; its Gaussian probability is at least `c'/t`. On this balance event `A_0,A_1=Θ(t)` and `O≤B`. Hence the normalized area measure obeys

`|μ(B_0)−1/2|≤C/t`, `μ((B_0∪B_1)^c)≤C/t`.       (4)

The unbiased field probability of (4) is at least `c t^(−2q−1)`. On that event the total mass is `Θ(t)`, so the exact BGS bias `M^κ` gives

`P_sphere((4)) ≥ c t^(−2q−1+κ) = c t^(−a−1)`.       (5)

Pick the fixed radius `r` sufficiently small that the independent two-curve SLE₈ block-order probability `q_SLE(r)>0` in `LEMMA_FIXED_CLUMP_SLE.md`. This probability depends on `r` but not on `t` or the field. By `LEMMA_TWO_BLOB_SUFFICIENCY.md`, (4) and block order imply `E_(C'/t)`. Equation (5) therefore yields

`P(E_(C'/t))≥c' t^(−a−1)`,

and hence (1)–(2). The use of a fixed radius avoids any uniform small-radius GMC estimate. The argument establishes only a one-sided bound for the BGS continuum event; it neither identifies the full small-ball exponent nor transfers it to finite meanders.

## Sources and dependencies

- [Borga–Gwynne–Sun, Definition 2.2 and Section 2.4](https://arxiv.org/html/2207.02319): sphere bias, area marks, and fixed-point SLE visits.
- [Wong, Proposition A.1 and Lemmas 3.7–3.8](https://arxiv.org/html/1902.04054): pinned-GMC regular variation and smooth-covariance transfer, as checked in `1342_alpha_tree_entropy/PROOF_PINNED_GMC_LOCAL_TAIL.md`.
- `PROOF_OFFBLOCK_SET_IDENTITY.md`, `LEMMA_TWO_BLOB_SUFFICIENCY.md`, `LEMMA_FIXED_CLUMP_SLE.md`: deterministic and SLE event inclusions, independently reviewed in `1342_alpha_tree_entropy/REFEREE_FIXED_CLUMP_SLE.md`.
