# Non-D-finiteness of the band [-1,2] meander series

Status: proved from the reviewed band exponent theorem; independently reviewed with verdict PASS in 1352_extensions_review/REFEREE_E1.md. All polynomial differential equations below are homogeneous and nonzero. Adding or deleting the constant term has no effect on the conclusion.

## 1. Statement and exact inputs

Let a_n count order-n meanders with the labels of Paper C, Section 6 in {-1,0,1,2}. Let m be the positive root of 4m^4-2m-1=0, c=(1+2m)/(3+4m), theta=arccos(c), alpha=1+pi/theta, and Lambda=(1+2m)(1+m)^2. Then F(z)=sum a_n z^n is not D-finite over C(z). Consequently it is not algebraic over C(z).

The project input is precisely BG1 (650, Section 3), which combines the upper theorem of 632 as repaired/reviewed in 639 with the lower theorem of 643 as repaired/reviewed in 650:

For every epsilon>0 there exist positive constants C_epsilon,c_epsilon such that for every n>=2,

c_epsilon Lambda^n n^(-alpha-epsilon) <= a_n <= C_epsilon Lambda^n n^(-alpha+epsilon).   (1)

We use no full asymptotic equivalent, ratio limit, unique dominant singularity, or sharp constant. The remaining inputs are a_n nonnegative integers and the arithmetic theorem specified in Section 4. The weak exponential bound a_n<=Cat_n^2<=16^n is also immediate from the definition.

## 2. Arithmetic of alpha, with an irreducibility certificate

Set u=2c. Inverting u=2(1+2m)/(3+4m) gives m=(2-3u)/(4u-4). Direct substitution gives

(4u-4)^4 (4m^4-2m-1)=4P(u),
P(u)=113u^4-312u^3+312u^2-128u+16.

The denominator does not vanish because u=1 would give 2+4m=3+4m. Thus P(u)=0. P is primitive. Modulo 3, its monic associate is f(X)=X^4-X-1. This polynomial has values 2,2,1 at 0,1,2, respectively. The three monic irreducible quadratics over F_3 are X^2+1, X^2+X+2, X^2+2X+2; the remainders of f on division by them are -X, 1-X, 1-X. Hence f has no linear or quadratic factor and is irreducible. Gauss's lemma shows P is irreducible over Q.

The monic minimal polynomial of u is P/113. It has nonintegral constant coefficient 16/113, so u is not an algebraic integer. If theta/pi were rational, exp(i theta) would be a root of unity, and u=exp(i theta)+exp(-i theta) would be an algebraic integer. This contradiction proves theta/pi irrational. If alpha were rational, theta/pi=1/(alpha-1) would be rational. Therefore alpha is irrational. Transcendence, although known here, is unnecessary.

## 3. An elementary radial lemma for logarithmic coefficient exponents

Lemma. Suppose b_n>=0 and, for a positive R and real alpha, the bounds
c_epsilon n^(-alpha-epsilon)<=b_n R^n<=C_epsilon n^(-alpha+epsilon)
hold eventually for every epsilon>0. Write B(z)=sum b_n z^n. For any integer k>=0 satisfying k+1>alpha,

lim_(t down to 0) log B^(k)(R(1-t))/log(1/t)=k+1-alpha.   (2)

Proof. Put beta=k-alpha and take 0<epsilon<beta+1. For n>=max(N_epsilon,2k), the falling factorial (n)_k is between positive constants times n^k, including k=0. Also (1-t)^(-k) stays bounded above and below near zero. The derivative is therefore bounded, up to fixed factors and an additional bounded polynomial, by the sums

sum n^(beta-epsilon)(1-t)^n and sum n^(beta+epsilon)(1-t)^n.

For each fixed p>-1, sum_(n>=N) n^p(1-t)^n is bounded above and below by positive constants times t^(-p-1). For the lower bound restrict to 1/t<=n<=2/t. For the upper bound use (1-t)^n<=exp(-tn), split n into intervals of length ceil(1/t), and sum the resulting exponentially decaying bounds; the first interval is bounded using sum_(n<=1/t) n^p=O(t^(-p-1)). Omitting finitely many terms does not change these estimates.

Thus the positive derivative lies between constant multiples of t^(-(k+1-alpha-epsilon)) and t^(-(k+1-alpha+epsilon)), for small enough t. The bounded polynomial can be absorbed since the lower exponent is positive. Divide logarithms by log(1/t), then send epsilon down to zero. This proves (2). The same bounds give radius exactly R by taking n-th roots. QED.

## 4. The arithmetic differential-equation input

Use the Andre-Chudnovsky-Katz theorem in the precise numbered formulation of Fischler and Rivoal, *On the values of G-functions*, Comment. Math. Helv. 89 (2014), 313-341, Theorem 6, Section 4.1, pp. 328-329, DOI 10.4171/CMH/321:
https://ems.press/content/serial-article-files/43350?nt=1

For a G-function, a minimal differential operator over its number field is Fuchsian and has rational local exponents. In particular, on a slit disk around a finite singularity R, every solution is a finite sum

sum_(q,j) (R-z)^q (log(R-z))^j H_(q,j)(R-z),   (3)

where q is rational, j is a nonnegative integer, and the H_(q,j) are holomorphic at zero. The cited theorem states both the regularity and this convergent local form. It attributes the result to Andre, Chudnovsky and Katz. We cite this combined numbered theorem rather than inventing separate numbers in their original monographs.

If F were D-finite over Q(z), it would be a G-function: coefficients are integers, so their common denominators are 1; conjugates are themselves; the exponential bound supplies the required archimedean bound. These check every arithmetic growth condition in Definition 1 of that paper. Its differential-equation condition is exactly the temporary D-finiteness assumption.

D-finiteness over C(z) descends to Q(z) for a rational-coefficient formal power series. Indeed, fix the order and coefficient-degree bounds of a proposed complex polynomial differential equation. Its coefficients are finitely many unknowns, constrained by infinitely many homogeneous linear equations over Q obtained by equating coefficients. Their row space has finite dimension and therefore a finite rational basis drawn from these equations. A nonzero complex kernel implies a nonzero rational kernel. Clearing denominators gives a Q[z] equation. Thus it is enough to exclude D-finiteness over Q(z).

## 5. Rational radial exponents and contradiction

Assume F is D-finite and put R=1/Lambda. Take any integer k with k+1>alpha. By (1) and the lemma, F^(k) tends to positive infinity radially and has the limit (2).

The local expansion (3) for F remains of the same finite form after k differentiations: powers shift by integers, logarithm degrees remain nonnegative integers, and analytic functions remain analytic. Put s=R-z>0 and choose a positive integer d clearing all denominators of its rational powers. In the variable w=s^(1/d), the derivative has a convergent Laurent series with polynomial coefficients in log w, with only finitely many negative powers of w. Combine all terms with the same power. Since the derivative is not zero, there is a least exponent l with nonzero polynomial coefficient P_l. Its highest logarithm term dominates the other terms at that exponent. All higher exponents are negligible (convergence of the finite analytic expansions ensures their sum is O(w times a fixed power of abs(log w)) relative to the minimal power, before any fixed polynomial division).

Therefore for some rational q=l/d, nonnegative integer j and nonzero constant A,

F^(k)(R-s)=A s^q (log s)^j (1+o(1)).

In particular its positive radial logarithmic growth divided by log(1/t), with s=Rt, is -q, a rational number. Equation (2) says it equals k+1-alpha. Hence alpha is rational, contradicting Section 2. This proves the theorem. Finally every algebraic univariate series in characteristic zero is D-finite, so F is not algebraic. QED.

## 6. Why several dominant singularities cause no problem

The proof only evaluates a nonnegative coefficient series and its derivative on the interval 0<z<R. It never transfers contributions from different singularities to coefficients. The Abelian sandwich in Section 3 already incorporates all coefficients. Contributions on other rays cannot cancel this positive radial divergence. At the single point R, Section 5 explicitly combines all local terms before choosing the first nonzero rational power. No coefficient asymptotic expansion or aperiodicity assumption is needed.
