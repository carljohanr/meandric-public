# Global near-block upper exponent alpha minus one

Status: proved; independent PASS reviews in
1341_alpha_tail/REFEREE_GLOBAL_NEARBLOCK_UPPER_A.md and
1342_alpha_tree_entropy/REFEREE_GLOBAL_NEARBLOCK_UPPER_A.md.

Let g=(17−sqrt(145))/3, a=4/g=(17+sqrt(145))/12, and α=a+1.
Let R_δ denote the actual rooted near-block event for two independent
whole-plane SLE8 curves on the unit-area BGS sphere. Set F(δ)=P(R_δ).
For every ζ>0 there is a finite C_ζ such that

    F(δ) <= C_ζ δ^(a−ζ),       0<δ<1.                    (1)

The intended improvement is from an unspecified positive upper exponent
to every exponent below a≈2.42013288. The predicted exponent α≈3.42013288
remains one power higher. This is a continuum theorem, not a discrete
meander asymptotic or a numerical growth-constant bound.

## Inputs and scopes

All raw masses refer to one fixed BGS outer-circle gauge. Write T for
total raw physical area and μ for its normalization.

1. The global area-mark rank lemma in
   1341/PROOF_ADAPTIVE_ORDINARY_CLUMP_MARK_RANK.md gives fixed ρ>0 such
   that C={some chordal ball of radius <=ρ has μ-mass >=1−ε} satisfies
   P(R_δ intersect C)<=(3ε+ε²)F(δ). Its independent review is
   1342/REFEREE_ADAPTIVE_ORDINARY_MARK_RANK.md. Fix ε=1/20, so
   λ=3ε+ε²=61/400<1/2. The precise positive choice is immaterial.
2. The fixed-frontier rank lemma in
   1341/PROOF_MACROSCOPIC_FIRST_FRONTIER_RANK.md gives fixed d0>0 with
   P(R_δ,diam Γ<d0)<=F(δ)/2. Reviewed in
   1343/REFEREE_MACROSCOPIC_FIRST_FRONTIER_RANK.md.
3. PROOF_LOW_TOTAL_MACROSCOPIC_FRONTIER.md gives, for every fixed v<1
   and P>0,
   P(R_δ,T<=δ^(−v),diam Γ>=d0)<=C_(v,P)δ^P.
   This is the new raw-total truncated GMC / multipoint SLE argument.
4. The reviewed high-total two-cell bound in Section 1 of
   1341/NEARBLOCK_DIFFUSE_THREE_ROUTES.md gives, for every small η>0,

       P(T>t,C^c) <= C_η t^(−a+η),       t>=1.           (2)

   Indeed absence of a dominant small ball extracts two separated fixed
   mesh cells carrying fixed fractions of T. Moment thresholds are q=a−1
   for inserted cells and a for ordinary cells. The sphere bias refunds
   κ=a−2. The worst of the three pair types is two inserted cells, with
   exponent 2q−κ=a. If either is ordinary the exponent is at least a+1.
   Constants depend on fixed ρ and ε; the finite mesh union adds no δ
   power. The same fixed-gauge finite-sphere-atlas estimate covers ∞.
   Reviews: 1342/REFEREE_TWO_CELL_HIGH_MASS.md and root's prior audit.

Input 4 does not need the near-block event or half-area balance.
Input 3 does need the near-block event, both independent curves, and
macroscopic Γ. No independent events are multiplied in the assembly.

## Absorption

Partition R_δ by small Γ, by C, and then by total mass. The classes need
not be disjoint; a union bound suffices. Inputs 1--4 with t=δ^(−v) give

    F(δ) <= (1/2+λ)F(δ) + C_(v,P)δ^P
                         + C_η δ^(v(a−η)).              (3)

Since 1/2−λ=139/400>0, absorb the first term. Given 0<ζ<a, choose
fixed η>0 and v<1 sufficiently close to 1 that v(a−η)>a−ζ, and choose
P>a. The resulting constants are finite and independent of δ. This
proves (1). For ζ>=a the statement follows by enlarging the constant.
No uniformity as v approaches 1 is needed or claimed.

The reviewed lower bound F(δ)>=c δ^α and (1) place any logarithmic
exponent, if it exists, in [a,α]. They do not establish existence of the
limit. Explicitly liminf log F(δ)/log δ>=a and limsup<=α.

## Remaining sharp gate

The same argument, retaining the ordinary-cell improvement instead of
discarding it, reduces the α−ζ bound to the high-total pin-pin class:
two distinct insertion neighborhoods carry fixed fractions of total raw
area while the curves nearly agree. The two cell totals need not each
be half of total mass. Consequently a Gaussian half-balance window
cannot be appended without a further argument about internal curve
ordering, or a decomposition with controlled balance errors. The
endpoint low-T estimate at v=1 is unnecessary for a logarithmic sharp
exponent: a sharp high-T pin-pin bound uniform in δ at fixed t would
combine with v approaching 1 in the same way.
