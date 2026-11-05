# Reciprocal head-edge fields in the paired quadratic comparison

2026-10-11. Complete sufficient argument written here; independent mathematical review is pending. It uses the standing tape root-return identity and the already reviewed finite-cap adjoint argument. There is no claim of a new accepted global endpoint or literature novelty.

## 1. A finite quadratic certificate

Let A be the row adjacency of any finite principal restriction of the physical tape machine. Let W,Z be positive vertex weights. Set

    D_i=sqrt(Z_i/W_i),   B=D A D^-1,   S=(B+B^T)/2,
    u_i=sqrt(W_i Z_i).

For every unordered edge in the support of S, choose real eta_ij with eta_ji=-eta_ij. Define t_ij=(u_j/u_i) exp(eta_ij), so t_ji=1/t_ij. If

    c_i = sum_j S_ij t_ij <= R                         (1)

for every vertex, then lambda_max(S)<=R. Indeed, for each unordered edge,

    2 v_i v_j <= t_ij v_i^2 + t_ji v_j^2,

since their difference is (sqrt(t_ij)v_i-v_j/sqrt(t_ij))^2. Summing gives

    R sum_i v_i^2 - v^T S v
       = sum_edges S_ij (sqrt(t_ij)v_i-v_j/sqrt(t_ij))^2
         +sum_i (R-c_i)v_i^2 >=0.

The tape has no one-step self-loops. If a general matrix has them, include S_ii in c_i and keep the same formula for off-diagonal edges. This is a finite sum-of-squares certificate.

A real nonnegative Perron vector v of B gives rho(A)=rho(B)=v^T Sv/(v^T v)<=R, including reducible caps. The standard finite-cap passage bounds every length-n root return by R^n and gives mu<=R^2. No infinite-space diagonal similarity or operator norm is invoked.

Substitution into (1) gives the useful row formula

    c_i = (1/2) sum_j A_ij [W_j/W_i] exp(eta_ij)
         +(1/2) sum_j A_ji [Z_j/Z_i] exp(eta_ij).        (2)

The same eta from the current vertex to its neighbor appears in both sums. In the incoming sum it is NOT the eta in the direction of the incoming forward edge. When eta=0, (2) is exactly the existing paired forward/adjoint family.

## 2. Only true head edges receive fields

Set eta=0 on every edge changing tape length. For a head shift, the tape word is unchanged. Pad the physical word by r+1 zeros on each side for evaluating weights and fields. Choose an integer q<=r+1 and a positive multiplier table E on binary cut words of length 2q.

For a right shift from physical head h to h+1, using positions in this same padded word, put

    exp(eta_R) = E[word[h-q+1 : h+q+1]].                (3)

For a left shift from h to h-1, put

    exp(eta_L) = 1/E[word[h-q : h+q]].                  (4)

The intervals use Python's exclusive upper endpoint. A right shift and its actual left inverse read the IDENTICAL cut word and reciprocal multipliers, so exp(eta_ij)exp(eta_ji)=1 exactly. This is a statement about a real unchanged word, not independently rewired exterior letters.

An edit and a later four-move undo are not identified as reciprocal one-step edges. Every edit edge simply receives multiplier one in both directions of S. This distinction is essential because the physical directed edit rules are not their own transposes.

## 3. Locality and the incoming comparison

Take the ordinary positive gram/window weight W=G_F U on the zero-padded tape, and Z=G_(Fcomp) U(comp(window)) on that same zero-padded tape. Equivalently Z evaluates the W formula after complementing the entire padded word, including its padding. No complemented tape is asserted to be a physical charge-one state.

The fields (3),(4) are visible in the standard extended row of length 2r+3 because q<=r+1. The ordinary gram locality identity and its incoming version give every ratio in (2). Actual cap or endpoint restrictions remove positive terms. Thus checking the complete formal incoming/outgoing row sums is sufficient for every actual finite restriction, including short tapes and genuine endpoints.

For incoming head shifts, eta is evaluated from the ORIGINAL physical row and the direction from the current head to the incoming predecessor. The complemented row is only a convenient way of evaluating the incoming Z weight ratio. Evaluating eta on a complemented or incoming-forward source without justification would be a different, potentially invalid certificate.

For this gate impose exact identities

    E(comp(cut))=E(cut),       E(reverse(cut))=1/E(cut).  (5)

Complement-even fields make the two formal row costs at z and comp(z) identical, so the usual paired grouping remains valid. Reflection-odd log fields preserve spatial reflection of the row problem. There are two free E parameters at q2 and twelve at q3; cuts fixed by a conflicting symmetry receive E=1. Signed lookup tables construct all other values as a free positive rational parameter or its exact reciprocal.

The row cost is a positive sum of monomials in F,U,E. Consequently its logarithm is convex in their logarithms. The tiny epigraph search proposes weights. Its final complete formal rows are rescored with Fraction arithmetic, preserving field reciprocity exactly.

## 4. This field is not a renamed vertex potential

The accepted q3 witness has a genuine nonunit field circulation on the legal root cycle

    (101,head1) --split100--> (11001,head2)
      --L--> (11001,head1) --merge110--> (101,head0)
      --R--> (101,head1).

Here the head indices are zero based. Its product of edge multipliers is exactly

    42295681667/34359738368 = 1.230966348288348... !=1.

Every ratio of a single positive vertex potential telescopes to one around this closed path. Therefore these edge multipliers cannot be represented as a vertex-potential ratio, at any radius. In particular the ordinary minimum-over-compatible-memory collapse does not directly apply to this parameterization.

This does not prove that no other vertex model can achieve the same numerical bound. The finite matrix always has a Perron supersolution; that supersolution need not share this local storage. The edge multipliers are quadratic-inequality parameters, not transition probabilities, so a nonunit circulation is allowed.

## 5. Relation to the lower-flow dual

Append the signed head-cut features to the existing log-monomial feature matrix. The same entropy-dual derivation as in 1001 adds zero total flux for each new signed head-cut feature. These are extra reciprocal-edge balance constraints on the surrogate dual flow.

The exact old r2 binding-flow witness violates these constraints: maximum absolute added flux is about0.004214 at q2 and0.024389 at q3. This identifies a defect that the enlarged comparison removes. It is a diagnostic of the finite dual relaxation, not a physical probability of an impossible path and not a lower-bound certificate.

## 6. Controls and limits

`verify.py` imports none of the producer. It uses direct literal incoming rules, checks every formal rational row, verifies all cut symmetry identities, and checks the actual cap rows plus true head-edge inverses at caps J2,J3,J4. The cap audit also includes the extra length-one state [1], strengthening the endpoint test. `curl.py` checks the exact legal cycle circulation above.

The tiny gates improve the exact r2/k2 paired optimum12.74928850..12.74928851 to q2 candidate12.739115792380378 and q3 candidate12.64013470314917. Both optimization runs terminate successfully and all final rational rows pass, but no optimum dual bracket for either enlarged family is claimed. Mathematical review and a fair optimized resource comparison are still required before scaling. The global upper record remains12.36524829705025.
