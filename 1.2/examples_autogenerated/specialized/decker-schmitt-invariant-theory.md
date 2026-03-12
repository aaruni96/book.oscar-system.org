---
layout: page
title: Invariant Theory
---

## Authors: Wolfram Decker, Lakshmi Ramesh and Johannes Schmitt

This page contains the code samples from the corresponding chapter in the OSCAR book. You can access the full chapter [here](https://link.springer.com/chapter/10.1007/978-3-031-62127-7_12).

### Invariants of symmetric groups

The invariants of the symmetric groups $S_4$ and $S_3$ serve as basic examples, see Examples 2 and 7.

```julia
julia> RS4 = invariant_ring(GF(2), symmetric_group(4))
Invariant ring
  of Sym(4)

julia> L = fundamental_invariants(RS4)
4-element Vector{MPolyDecRingElem{FqFieldElem, FqMPolyRingElem}}:
 x[1] + x[2] + x[3] + x[4]
 x[1]*x[2] + x[1]*x[3] + x[2]*x[3] + x[1]*x[4] + x[2]*x[4] + x[3]*x[4]
 x[1]*x[2]*x[3] + x[1]*x[2]*x[4] + x[1]*x[3]*x[4] + x[2]*x[3]*x[4]
 x[1]*x[2]*x[3]*x[4]

julia> is_algebraically_independent(L)
true

julia> RS3 = invariant_ring(QQ, symmetric_group(3));

julia> R = polynomial_ring(RS3)
Multivariate polynomial ring in 3 variables over QQ graded by
  x[1] -> [1]
  x[2] -> [1]
  x[3] -> [1]

julia> x = gens(R);

julia> reynolds_operator(RS3, x[1]^2)
1//3*x[1]^2 + 1//3*x[2]^2 + 1//3*x[3]^2
```

### The Heisenberg group of level 3

The Heisenberg group of level 3 in its Schrödinger representation is studied in Examples 18, 21, and 23.

```julia
julia> K, a = cyclotomic_field(3, "a");

julia> S = matrix(K, [0 0 1; 1 0 0; 0 1 0])
[0   0   1]
[1   0   0]
[0   1   0]

julia> T = matrix(K, [1 0 0; 0 a 0; 0 0 -a-1])
[1   0        0]
[0   a        0]
[0   0   -a - 1]

julia> H3 = matrix_group(S, T)
Matrix group of degree 3
  over cyclotomic field of order 3

julia> order(H3)
27

julia> exponent(H3)
3

julia> is_abelian(H3)
false

julia> R, (x, y, z) = graded_polynomial_ring(K, ["x", "y", "z"]);

julia> RH3 = invariant_ring(R, H3);

julia> MSH3 = molien_series(RH3)
(-t^6 + t^3 - 1)//(t^9 - 3*t^6 + 3*t^3 - 1)

julia> expand(MSH3, 10)
1 + 2*t^3 + 4*t^6 + 7*t^9 + O(t^11)

julia> P = basis(RH3, 3)
2-element Vector{MPolyDecRingElem{AbsSimpleNumFieldElem, AbstractAlgebra.Generic.MPoly{AbsSimpleNumFieldElem}}}:
 x*y*z
 x^3 + y^3 + z^3

julia> MP = minimal_primes(ideal(R, P));

julia> length(MP)
9

julia> [[ kernel(reduce(hcat, [ K[coeff(p, x); coeff(p, y); coeff(p, z)] for p in gens(Q) ])) for Q in MP ]]
1-element Vector{Vector{AbstractAlgebra.Generic.MatSpaceElem{AbsSimpleNumFieldElem}}}:
 [[a+1 0 1], [-a 0 1], [-1 0 1], [0 a+1 1], [0 -a 1], [0 -1 1], [-a 1 0], [a+1 1 0], [-1 1 0]]

julia> S, s = QQ["t"]; T = fraction_field(S); t = T(s);

julia> RR, (X, Y, Z) = graded_polynomial_ring(T, [ "X", "Y", "Z"]);

julia> F = X^3+Y^3+Z^3;

julia> G = t*F+hessian(F)
t*X^3 + 216*X*Y*Z + t*Y^3 + t*Z^3

julia> L = syzygy_generators([hessian(G), F, hessian(F)]);

julia> collect(coefficients(L[1]))
3-element Vector{AbstractAlgebra.Generic.FracFieldElem{QQPolyRingElem}}:
 -1
 -279936*t
 t^3 + 93312

julia> C, iC = center(H3);

julia> gens(C)[1]
[-a - 1        0        0]
[     0   -a - 1        0]
[     0        0   -a - 1]

julia> Q, pQ = quo(H3, C)
(Pc group of order 9, Hom: H3 -> Q)

julia> describe(Q)
"C3 x C3"

julia> SSG = filter(cls -> order(representative(cls)) == 3 && length(cls) == 3, subgroup_classes(H3));

julia> length(SSG)
4

julia> L1 = [basis(invariant_ring(R, H), 1)[1] for H in collect(SSG[1])]
3-element Vector{MPolyDecRingElem{AbsSimpleNumFieldElem, AbstractAlgebra.Generic.MPoly{AbsSimpleNumFieldElem}}}:
 x + y + z
 x + (-a - 1)*y + a*z
 x + a*y + (-a - 1)*z

julia> Triangle1 = L1[1]*L1[2]*L1[3]
x^3 - 3*x*y*z + y^3 + z^3
```

### Non-trivial module syzygies

This example of an invariant ring which is not Cohen--Macaulay is presented in Examples 19 and 26.

```julia
julia> G = permutation_group(4, [ cperm([1, 2, 3, 4]) ]);

julia> RG = invariant_ring(GF(2), G);

julia> primary_invariants(RG)
4-element Vector{MPolyDecRingElem{FqFieldElem, FqMPolyRingElem}}:
 x[1] + x[2] + x[3] + x[4]
 x[1]*x[3] + x[2]*x[4]
 x[1]*x[2] + x[2]*x[3] + x[1]*x[4] + x[3]*x[4]
 x[1]*x[2]*x[3]*x[4]

julia> secondary_invariants(RG)
5-element Vector{MPolyDecRingElem{FqFieldElem, FqMPolyRingElem}}:
 1
 x[1]^2*x[2] + x[2]^2*x[3] + x[3]^2*x[4] + x[1]*x[4]^2
 x[1]*x[2]*x[3] + x[1]*x[2]*x[4] + x[1]*x[3]*x[4] + x[2]*x[3]*x[4]
 x[1]^2*x[2]*x[3] + x[2]^2*x[3]*x[4] + x[1]*x[3]^2*x[4] + x[1]*x[2]*x[4]^2
 x[1]^3*x[2]*x[3] + x[2]^3*x[3]*x[4] + x[1]*x[3]^3*x[4] + x[1]*x[2]*x[4]^3

julia> M, MtoR, StoR = module_syzygies(RG);

julia> M
Subquotient of submodule with 5 generators
  1: e[1]
  2: e[2]
  3: e[3]
  4: e[4]
  5: e[5]
by submodule with 1 generator
  1: t2*e[2] + (t2 + t3)*e[3] + t1*e[4]
```

### The Klein four-group

The invariants of the Klein four-group in a regular representation are computed in Examples 22, 24, and 25 both in characteristic 0 and characteristic 2.

```julia
julia> A = matrix(QQ, [0 1 0 0; 1 0 0 0; 0 0 0 1; 0 0 1 0])
[0   1   0   0]
[1   0   0   0]
[0   0   0   1]
[0   0   1   0]

julia> B = matrix(QQ, [0 0 1 0; 0 0 0 1; 1 0 0 0; 0 1 0 0])
[0   0   1   0]
[0   0   0   1]
[1   0   0   0]
[0   1   0   0]

julia> K4 = matrix_group(A, B)
Matrix group of degree 4
  over rational field

julia> small_group_identification(K4)
(4, 2)

julia> describe(K4)
"C2 x C2"

julia> RK4 = invariant_ring(K4);

julia> MSK4 = molien_series(RK4)
(t^2 - t + 1)//(t^6 - 2*t^5 - t^4 + 4*t^3 - t^2 - 2*t + 1)

julia> expand(MSK4, 4)
1 + t + 4*t^2 + 5*t^3 + 11*t^4 + O(t^5)

julia> primary_invariants(RK4)
4-element Vector{MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}}:
 x[1] + x[2] + x[3] + x[4]
 x[1]^2 + x[2]^2 + x[3]^2 + x[4]^2
 x[1]*x[2] + x[3]*x[4]
 x[1]*x[3] + x[2]*x[4]

julia> secondary_invariants(RK4)
2-element Vector{MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}}:
 1
 x[1]^3 + x[2]^3 + x[3]^3 + x[4]^3

julia> fundamental_invariants(RK4)
5-element Vector{MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}}:
 x[1] + x[2] + x[3] + x[4]
 x[1]^2 + x[2]^2 + x[3]^2 + x[4]^2
 x[1]*x[2] + x[3]*x[4]
 x[1]*x[3] + x[2]*x[4]
 x[1]^3 + x[2]^3 + x[3]^3 + x[4]^3

julia> A = matrix(GF(2), [0 1 0 0; 1 0 0 0; 0 0 0 1; 0 0 1 0]);

julia> B = matrix(GF(2), [0 0 1 0; 0 0 0 1; 1 0 0 0; 0 1 0 0]);

julia> RK42 = invariant_ring(matrix_group(A, B));

julia> primary_invariants(RK42)
4-element Vector{MPolyDecRingElem{FqFieldElem, FqMPolyRingElem}}:
 x[1] + x[2] + x[3] + x[4]
 x[2]*x[3] + x[1]*x[4]
 x[1]*x[2] + x[1]*x[3] + x[2]*x[4] + x[3]*x[4]
 x[1]*x[2]*x[3]*x[4]

julia> secondary_invariants(RK42)
4-element Vector{MPolyDecRingElem{FqFieldElem, FqMPolyRingElem}}:
 1
 x[1]*x[3] + x[2]*x[4]
 x[2]^2*x[3] + x[2]*x[3]^2 + x[1]^2*x[4] + x[1]*x[4]^2
 x[1]*x[2]^2*x[3]^2 + x[1]*x[2]*x[3]^3 + x[1]^3*x[3]*x[4] + x[2]^3*x[3]*x[4] + x[2]^2*x[3]^2*x[4] + x[1]^2*x[2]*x[4]^2 + x[1]^2*x[3]*x[4]^2 + x[1]*x[2]*x[4]^3

julia> fundamental_invariants(RK42)
6-element Vector{MPolyDecRingElem{FqFieldElem, FqMPolyRingElem}}:
 x[1] + x[2] + x[3] + x[4]
 x[1]*x[3] + x[2]*x[4]
 x[2]*x[3] + x[1]*x[4]
 x[1]*x[2] + x[1]*x[3] + x[2]*x[4] + x[3]*x[4]
 x[2]^2*x[3] + x[2]*x[3]^2 + x[1]^2*x[4] + x[1]*x[4]^2
 x[1]*x[2]*x[3]*x[4]
```

### A complex reflection group

The fundamental invariants of a complex reflection group which is connected to an application in coding theory are discussed in Example 27.

```julia
julia> Qt, t = QQ["t"];

julia> K, (a, b) = number_field([t^2 - 2, t^2 + 1], ["a", "b"]);

julia> M1 = 1/a*matrix(K, [1 1; 1 -1])
[1//2*a    1//2*a]
[1//2*a   -1//2*a]

julia> M2 = matrix(K, [1 0; 0 b])
[1   0]
[0   b]

julia> G = matrix_group(M1, M2);

julia> order(G)
192

julia> RG = invariant_ring(G);

julia> f = fundamental_invariants(RG)
2-element Vector{MPolyDecRingElem{AbsNonSimpleNumFieldElem, AbstractAlgebra.Generic.MPoly{AbsNonSimpleNumFieldElem}}}:
 x[1]^8 + 14*x[1]^4*x[2]^4 + x[2]^8
 x[1]^20*x[2]^4 - 4*x[1]^16*x[2]^8 + 6*x[1]^12*x[2]^12 - 4*x[1]^8*x[2]^16 + x[1]^4*x[2]^20

julia> A, AtoR = affine_algebra(RG)
(Quotient of multivariate polynomial ring by ideal (0), Hom: A -> graded multivariate polynomial ring)

julia> modulus(A)
Ideal generated by
  0
```

### A pseudo-reflection group in positive characteristic

An example of a pseudo-reflection group whose invariant ring is not a polynomial ring is presented in Example 28.

```julia
julia> K = GF(3);

julia> M1 = matrix(K, [1 0 1 0; 0 1 0 0; 0 0 1 0; 0 0 0 1])
[1   0   1   0]
[0   1   0   0]
[0   0   1   0]
[0   0   0   1]

julia> M2 = matrix(K, [1 0 0 0; 0 1 0 1; 0 0 1 0; 0 0 0 1])
[1   0   0   0]
[0   1   0   1]
[0   0   1   0]
[0   0   0   1]

julia> M3 = matrix(K, [1 0 1 1; 0 1 1 1; 0 0 1 0; 0 0 0 1])
[1   0   1   1]
[0   1   1   1]
[0   0   1   0]
[0   0   0   1]

julia> G = matrix_group(M1, M2, M3);

julia> all(is_pseudo_reflection, gens(G))
true

julia> RG = invariant_ring(G);

julia> fundamental_invariants(RG)
5-element Vector{MPolyDecRingElem{FqFieldElem, FqMPolyRingElem}}:
 x[3]
 x[4]
 x[1]^3*x[3] + 2*x[1]*x[3]^3 + x[2]^3*x[4] + 2*x[2]*x[4]^3
 x[1]^9 + 2*x[1]^3*x[3]^6 + 2*x[1]^3*x[3]^4*x[4]^2 + x[1]*x[3]^6*x[4]^2 + 2*x[1]^3*x[3]^2*x[4]^4 + x[1]*x[3]^4*x[4]^4 + 2*x[1]^3*x[4]^6 + x[1]*x[3]^2*x[4]^6
 x[2]^9 + 2*x[2]^3*x[3]^6 + x[1]^3*x[3]^5*x[4] + 2*x[1]*x[3]^7*x[4] + x[2]*x[3]^6*x[4]^2 + x[1]^3*x[3]^3*x[4]^3 + 2*x[1]*x[3]^5*x[4]^3 + 2*x[2]^3*x[4]^6

julia> A, AtoR = affine_algebra(RG)
(Quotient of multivariate polynomial ring by ideal (y1^6*y2^2*y3 + y1^3*y4 + y2^3*y5 + 2*y3^3), Hom: A -> graded multivariate polynomial ring)

julia> modulus(A)
Ideal generated by
  y1^6*y2^2*y3 + y1^3*y4 + y2^3*y5 + 2*y3^3
```

### Invariants of linearly reductive groups

Computations with invariants of linearly reductive groups can be found in Examples 34, 35, and 36.

```julia
julia> G = linearly_reductive_group(:SL, 2, QQ)
Reductive group SL2
  over QQ

julia> r = representation_on_forms(G, 2)
Representation of SL2
  on symmetric forms of degree 2

julia> RG = invariant_ring(r)
Invariant Ring of
graded multivariate polynomial ring in 3 variables over QQ
  under group action of SL2

julia> fundamental_invariants(RG)
1-element Vector{MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}}:
 -X[1]*X[3] + X[2]^2

julia> G = linearly_reductive_group(:SL, 3, QQ);

julia> r = representation_on_forms(G, 3);

julia> S, x = graded_polynomial_ring(QQ, "x" => 1:10);

julia> RG = invariant_ring(S, r);

julia> 75*reynolds_operator(RG, x[5]^4)
x[1]*x[4]*x[8]*x[10] - x[1]*x[4]*x[9]^2 - x[1]*x[5]*x[7]*x[10] + x[1]*x[5]*x[8]*x[9] + x[1]*x[6]*x[7]*x[9] - x[1]*x[6]*x[8]^2 - x[2]^2*x[8]*x[10] + x[2]^2*x[9]^2 + x[2]*x[3]*x[7]*x[10] - x[2]*x[3]*x[8]*x[9] + x[2]*x[4]*x[5]*x[10] - x[2]*x[4]*x[6]*x[9] - 2*x[2]*x[5]^2*x[9] + 3*x[2]*x[5]*x[6]*x[8] - x[2]*x[6]^2*x[7] - x[3]^2*x[7]*x[9] + x[3]^2*x[8]^2 - x[3]*x[4]^2*x[10] + 3*x[3]*x[4]*x[5]*x[9] - x[3]*x[4]*x[6]*x[8] - 2*x[3]*x[5]^2*x[8] + x[3]*x[5]*x[6]*x[7] + x[4]^2*x[6]^2 - 2*x[4]*x[5]^2*x[6] + x[5]^4

julia> T = torus_group(QQ,2)
Torus of rank 2
  over QQ

julia> r = representation_from_weights(T, [1 0; 0 1; -1 -1; -1 1]);

julia> RT = invariant_ring(r);

julia> fundamental_invariants(RT)
2-element Vector{MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}}:
 X[1]*X[2]*X[3]
 X[1]^2*X[3]*X[4]
```

### The Cox ring of a linear quotient

This computation of a Cox ring of a linear quotient is worked out in Example 40.

```julia
julia> K, z_3 = cyclotomic_field(3)
(Cyclotomic field of order 3, z_3)

julia> g1 = matrix(K, 3, 3, [ -1 0 0; 0 1 0; 0 0 1 ]);

julia> g2 = matrix(K, 3, 3, [ 0 0 1; 1 0 0; 0 1 0 ]);

julia> G = matrix_group(g1, g2)
Matrix group of degree 3
  over cyclotomic field of order 3

julia> chi = character_table(G)[6];

julia> RG = invariant_ring(G)
Invariant ring
  of matrix group of degree 3 over K

julia> semi_invariants(RG, chi)
2-element Vector{MPolyDecRingElem{AbsSimpleNumFieldElem, AbstractAlgebra.Generic.MPoly{AbsSimpleNumFieldElem}}}:
 x[1]^2 + z_3*x[2]^2 + (-z_3 - 1)*x[3]^2
 x[1]^2*x[2]^2 + (-z_3 - 1)*x[1]^2*x[3]^2 + z_3*x[2]^2*x[3]^2

julia> H, HtoG = subgroup_of_pseudo_reflections(G)
(Matrix group of degree 3 over K, Hom: H -> G)

julia> RH = invariant_ring(H)
Invariant ring
  of matrix group of degree 3 over K

julia> A, AtoR = affine_algebra(RH)
(Quotient of multivariate polynomial ring by ideal (0), Hom: A -> graded multivariate polynomial ring)

julia> map(AtoR, gens(A))
3-element Vector{MPolyDecRingElem{AbsSimpleNumFieldElem, AbstractAlgebra.Generic.MPoly{AbsSimpleNumFieldElem}}}:
 x[1]^2
 x[2]^2
 x[3]^2

julia> VG = linear_quotient(G)
Linear quotient by matrix group of degree 3 over K

julia> B, BtoR = cox_ring(VG)
(Quotient of multivariate polynomial ring by ideal (), Hom: B -> graded multivariate polynomial ring)

julia> map(BtoR, gens(B))
3-element Vector{MPolyDecRingElem{AbsSimpleNumFieldElem, AbstractAlgebra.Generic.MPoly{AbsSimpleNumFieldElem}}}:
 z_3*x[1]^2 + (-z_3 - 1)*x[2]^2 + x[3]^2
 x[1]^2 + x[2]^2 + x[3]^2
 (-z_3 - 1)*x[1]^2 + z_3*x[2]^2 + x[3]^2

julia> f1 = BtoR(B[1]); f2 = BtoR(B[2]); f3 = BtoR(B[3]);

julia> f1^G(g2) == z_3*f1
true

julia> f2^G(g2) == f2
true

julia> f3^G(g2) == z_3^(-1)*f3
true

julia> grading_group(B)
Z/3

julia> degree(B[1])
Abelian group element [2]

julia> degree(B[2])
Abelian group element [0]

julia> degree(B[3])
Abelian group element [1]

julia> is_zero(modulus(B))
true
```
