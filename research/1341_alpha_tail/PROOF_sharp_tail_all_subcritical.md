# Sharp APPS cell-count exponent for every c<1, assembled from barrier S1 and S2

**Status.** Proved and reviewed. The sharp upper, deterministic-square descendant lemma, and general-c<1 upper-barrier S1 were independently reviewed in `1343_alpha_meander_bridge/REFEREE_SHARP_UPPER.md`, `1342_alpha_tree_entropy/REFEREE_APPS_DESCENDANTS.md`, and `1343_alpha_meander_bridge/REFEREE_SPINE_GENERAL.md`. Root audited the marked-count assembly in `1340_alpha_campaign/REFEREE_SHARP_TAIL_ASSEMBLY.md`. This supersedes the restricted `PROOF_sharp_tail_c_le_minus2.md`.

**Theorem.** Fix ε>0 and Q>2 (equivalently c<1). Put γ=Q−√(Q²−4). For the APPS terminal-cell count N,

    lim_(n→∞) log P_Q(N≥n)/log n = −4/γ².

**Proof.** The reviewed upper proof gives the limsup ≤−4/γ². For the lower bound set L=log 2, d=√(Q²−4), a=Q+d, and choose depth j→∞. Put y_j=djL−Q log ε. Fix δ>0. For every deterministic depth-j square S define B_S={h_S∈[ajL,ajL+1]} and the ceiling C_S={h_{S_l}≤alL+δjL for every ancestor at depth l≤j}. Let Z_j count all B_S squares, T_j count B_S∩C_S squares, and G_j count B_S∩C_S squares that are reachable and active and have at least exp((γ−η)y_j) active descendants along their whole paths, for an arbitrary η>0.

The uniform ancestral-survival estimate in `PROOF_spine_hit.md` gives P(S reachable and active | B_S)≥η₀>0, and its conditional ceiling-failure bound is o(1). The descendant lemma in `1340_alpha_campaign/PROOF_DESCENDANTS.md` gives, uniformly in S,

    P(S has ≥exp((γ−η)y_j) active-path descendants | B_S)=1−o(1).

These statements need no conditional independence. The union bound gives P(C_S, S reachable and descendant-successful | B_S)≥η₀−o(1). Consequently E G_j≥(η₀−o(1))E Z_j. Since G_j≤T_j, the barrier S1 second-moment estimate yields

    P(G_j>0)≥(E G_j)²/E G_j²
             ≥(η₀−o(1))²(E Z_j)²/E T_j²
             ≥exp(−(a+aδ/d+o(1))y_j).

On G_j>0, the subdivision has at least exp((γ−η)y_j) distinct active squares at one depth, each with a terminal descendant, so N≥exp((γ−η)y_j). Monotonicity extends the resulting bound from these exponentially spaced thresholds to all n; the bounded spacing in log n changes only O(1). Thus

    liminf_(n→∞) log P(N≥n)/log n ≥−(a+aδ/d)/(γ−η).

Let δ↓0 and η↓0. As aγ=(Q+d)(Q−d)=4, a/γ=4/γ². This matches the limsup. ∎

**Review boundary.** The reviews cover the upper-barrier second moment, the descendant lemma's uniform negative-moment input, and the marked-count assembly. The marking step does not condition on a first-found square or on an exploratory σ-field.

## Compensated cumulative determinant moment at base c=0

The theorem above and the exact Gaussian tilt in `PROOF_tilt.md` give a determinant statement. At base central charge c=0, set Q₀²=25/6 and for real k>−1/2 set Q_k²=Q₀²+k/3, so c_k=−2k<1. Write `D=exp(−Σ_S x_S²/6)` on `{N=m}`, and let `A₀=det′_ζΔ_g0/Vol_g0`. The Polyakov–Alvarez formula gives the PA-defined determinant `det′_PA Δ_g=A₀m^(1/3)D`. Define the compensated cumulative moment

    C_k(n)=Σ_(m≥n)(Q_k/Q₀)^m A₀^(−k)m^(−k/3)
                         E_(Q₀)[(det′_PA Δ_g)^k;N=m].

Then the exact finite-dimensional Gaussian tilt identity gives, term by term,

    C_k(n)=P_(Q_k)(N≥n).

Hence, for every real k>−1/2,

    lim_(n→∞) log C_k(n)/log n=−x_k,
    x_k=(13+2k+sqrt((13+2k)²−144))/12.

In particular, k=2 yields `x₂=(17+sqrt145)/12`, the c=−4 APPS exponent. The endpoint k=−1/2 corresponds to Q=2 and is not covered. Here `det′_PA Δ_g` is the APPS determinant **defined by the Polyakov–Alvarez formula** for its distributional conformal metric; `A₀=det′_ζΔ_g0/Vol_g0` uses the ordinary smooth reference metric. No claim about spectral-zeta continuation on the nonsmooth metric is made. The polynomial `m^(−k/3)` compensation is essential for the PA-defined determinant; without it the displayed exact identity does not hold. This is a continuum APPS theorem only. It neither proves the meander determinant limit nor a local coefficient asymptotic.
