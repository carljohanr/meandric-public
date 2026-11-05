# Explicit global five-sevenths gain from the averaged-target bootstrap

Status: proved; independent PASS in
1343/REFEREE_GLOBAL_FIVE_SEVENTHS_GAIN.md. The target degree-ten appendix
has independent review in 1341/REFEREE_GLOBAL_AVERAGED_ALL_COLOR_POWER.md.

Let a=(17+sqrt(145))/12 and q=a-1. The conclusion is

    F(delta) <= C delta^(a+5/7),

for the same rooted BGS continuum half-prefix event. The exponent is
approximately 3.13441860. No discrete assertion is made.

## Named inputs

1. 1343/CANDIDATE_GLOBAL_AVERAGED_ALL_COLOR_POWER.md gives the SAME completed
   sphere Good functional an averaged bound epsilon^(a+H)d^(-10), for
   H=7/10, with its explicit degree-ten audit independently passed. Its all-prefix power
   and finite-degree claim already have independent 1341/1342 reviews.
2. The older reviewed averaged bound supplies epsilon^a d^(-q), up to
   arbitrarily small strict losses.
3. LEMMA_HYBRID_TARGET_RANK_BOOTSTRAP.md, independently reviewed in
   1342/REFEREE_HYBRID_TARGET_RANK_BOOTSTRAP.md, transfers these two target
   bounds to a stronger retained local rank-integral power.
4. The reviewed all-color first-light coarea and global partition in
   CANDIDATE_ALL_COLOR_LOCAL_TRACE_SLICING.md transfer any such strict
   retained power below one. Their use does not require that this power
   be below the old a/(a+1); that value came solely from the old padding
   calculation, now replaced by input 3.

## Exact parameter margins

Use A=a+7/10, D=10, xi=1/40. The limiting old-small-rank exponent is

    lambda_old=(a+1/40)/(a+1).

It exceeds 5/7 exactly when a>193/80. Since sqrt(145)>12,
a>29/12>193/80. The limiting new-large-rank exponent is

    lambda_new=[A-(10-A)/40]/(A+1).

It exceeds 5/7 exactly when A>90/29. This follows from
A>29/12+7/10=187/60>90/29. The new-cutoff condition
xi<A/(D+1) holds because A>3, hence A/11>3/11>1/40.
Discarded ranks may use L0=3/4>5/7. Also H=7/10 is permitted by
the averaged theorem: a>29/12>7/3 implies a/(a+1)>7/10.

Choose a strict lambda between 5/7 and min(lambda_old,lambda_new,L0).
All CM, padding-moment, target, and curve-coefficient losses can be
chosen sufficiently small after this fixed positive gap. The hybrid
lemma gives the absolute one-pin estimate

 E[1_{B>=t} L_D(delta)]
 <= C_eta t^(-q+eta)[delta^(1-eta)+t^(-lambda)].

The weak disk is chosen after these finite Holder parameters, before
any t,delta or mesh choices. No endpoint limit of constants is used.

## Transfer and exhaustive global sum

In the all-color first-light proof, replace WG7 by this new one-pin
bound. Conditional strong-field integration and Holder still permit
p arbitrarily close to q and the Holder exponent arbitrarily close to
one. The same coefficient rate r>sigma>5/8 makes every first-light
scale summable. At the initial slice, retain the strong-core p-moment
and apply the same Holder argument without an m factor; this avoids
invoking an old two-pin statement with its former lambda restriction.
The strong moment coefficient has every fixed finite exterior moment
because the entire fixed strong Dirichlet field has already been
integrated. Thus both slice zero and subsequent slices have protected
high-band exponents exceeding a+5/7 after strict losses.

The weak retained trace law and both root orientations are unchanged.
Remote-tip source classes and ordinary-cell exceptions have every
strict gain below one. Choose common high-band threshold exponent v<1
sufficiently close to one after the finite geometry and all strict
margins. The same exhaustive pin/tip partition and low-total theorem
then yield

    (139/400)F(delta) <= C delta^(a+5/7).

This proves the conclusion. Both the degree-ten target input and this
substitution have passed independent review. The sharp a+1 target remains open.
