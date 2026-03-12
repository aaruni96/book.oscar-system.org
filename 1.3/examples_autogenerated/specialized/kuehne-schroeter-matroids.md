---
layout: page
title: Matroids
---

## Authors: Daniel Corey, Lukas Kühne and Benjamin Schröter

This page contains the code samples from the corresponding chapter in the OSCAR book. You can access the full chapter [here](https://link.springer.com/chapter/10.1007/978-3-031-62127-7_14).

Interested readers can find more information in [this article](https://arxiv.org/abs/2311.08792).

```julia
julia> B = [[1,2],[1,3],[1,4],[2,3],[2,4]];

julia> M = matroid_from_bases(B,4)
Matroid of rank 2 on 4 elements

julia> independent_sets(M)
10-element Vector{Vector{Int64}}:
 []
 [1]
 [2]
 [3]
 [4]
 [1, 4]
 [2, 4]
 [1, 3]
 [2, 3]
 [1, 2]

julia> flats(M)
5-element Vector{Vector{Int64}}:
 []
 [1]
 [2]
 [3, 4]
 [1, 2, 3, 4]

julia> circuits(M)
3-element Vector{Vector{Int64}}:
 [3, 4]
 [1, 2, 4]
 [1, 2, 3]

julia> X = matrix(GF(2), [1 0 1 0 1 0 1; 0 1 1 0 0 1 1; 0 0 0 1 1 1 1])
[1   0   1   0   1   0   1]
[0   1   1   0   0   1   1]
[0   0   0   1   1   1   1]

julia> F = matroid_from_matrix_columns(X)
Matroid of rank 3 on 7 elements

julia> is_isomorphic(F, fano_matroid())
true

julia> G = automorphism_group(F);

julia> describe(G)
"PSL(3,2)"
```

```julia
julia> nonbases(vamos_matroid())
5-element Vector{Vector{Int64}}:
 [1, 2, 3, 4]
 [1, 2, 5, 6]
 [1, 2, 7, 8]
 [3, 4, 5, 6]
 [5, 6, 7, 8]

julia> M = moebius_kantor_matroid();

julia> nonbases(M)
8-element Vector{Vector{Int64}}:
 [1, 2, 3]
 [1, 4, 5]
 [1, 7, 8]
 [2, 4, 6]
 [2, 5, 8]
 [3, 4, 7]
 [3, 6, 8]
 [5, 6, 7]

julia> realization_space(M)
The realization space is
  [1   0   1   0   1     0         1    1]
  [0   1   1   0   0     1         1   x1]
  [0   0   0   1   1   -x1   -x1 + 1    1]
in the multivariate polynomial ring in 1 variable over ZZ
within the vanishing set of the ideal
Ideal (x1^2 - x1 + 1)
avoiding the zero loci of the polynomials
RingElem[x1, x1 - 1, x1^2 + 1]

julia> X = matrix(GF(2), [1 0 1 0 1 0 1; 0 1 1 0 0 1 1; 0 0 0 1 1 1 1])
[1   0   1   0   1   0   1]
[0   1   1   0   0   1   1]
[0   0   0   1   1   1   1]

julia> F = matroid_from_matrix_columns(X)
Matroid of rank 3 on 7 elements

julia> realization_space(F, simplify=false)
The realization space is
  [0   1    1    1    1   0   0]
  [1   0    1   x2    0   1   0]
  [1   0   x1    0   x3   0   1]
in the multivariate polynomial ring in 3 variables over ZZ
within the vanishing set of the ideal
Ideal (2, x3 + 1, x2 - 1, x1 + 1)
avoiding the zero loci of the polynomials
RingElem[x2, x3, x1, x1*x2 - x2*x3 + x3, x1 - x3 - 1, x1 + x2 - 1]

julia> realization_space(F, simplify=true)
The realization space is
  [0   1   1   1   1   0   0]
  [1   0   1   1   0   1   0]
  [1   0   1   0   1   0   1]
in the integer ring
within the vanishing set of the ideal
2ZZ

julia> is_realizable(F, char=5)
false

julia> realization_space(non_fano_matroid())
The realization space is
  [1   1   0   0   1   1   0]
  [0   1   1   1   1   0   0]
  [0   1   1   0   0   1   1]
in the integer ring
avoiding the zero loci of the polynomials
RingElem[2]

julia> is_realizable(vamos_matroid())
false

julia> is_realizable(M, q=9)
true

julia> qs = [2,3,4,5,7,8,9,11,13];

julia> println([is_realizable(M, q=a) for a in qs])
Bool[0, 1, 1, 0, 1, 0, 1, 0, 1]

julia> M = cycle_matroid(complete_graph(4));

julia> realization_space(M)
The realization space is
  [1   0   1   0   1    0]
  [0   1   1   0   0    1]
  [0   0   0   1   1   -1]
in the integer ring

julia> db = Polymake.Polydb.get_db();

julia> collection = db["Matroids.Small"];

julia> query = Dict("RANK"=>3,"N_ELEMENTS"=>9,"SIMPLE"=>true);

julia> results = Polymake.Polydb.find(collection, query);

julia> oscar_matroids = [Matroid(pm) for pm in results];

julia> length(oscar_matroids)
383

julia> char_0_matroids = [M for M in oscar_matroids 
        if is_realizable(M, char=0)];

julia> length(char_0_matroids)
370

julia> RS = realization_space(pappus_matroid(), char=0)
The realization space is
  [1   0   1   0   x2   x2                 x2^2    1    0]
  [0   1   1   0    1    1   -x1*x2 + x1 + x2^2    1    1]
  [0   0   0   1   x2   x1                x1*x2   x1   x2]
in the multivariate polynomial ring in 2 variables over QQ
avoiding the zero loci of the polynomials
RingElem[x1 - x2, x2, x1, x2 - 1, x1 + x2^2 - x2, x1 - 1, x1*x2 - x1 - x2^2]

julia> realization(RS)
One realization is given by
  [1   0   1   0   2   2   4   1   0]
  [0   1   1   0   1   1   1   1   1]
  [0   0   0   1   2   3   6   3   2]
in the rational field

julia> results = collection = db = 0; # allow connection cleanup
```

```julia
julia> M = cycle_matroid(complete_graph(4))
Matroid of rank 3 on 6 elements

julia> tutte_polynomial(M)
x^3 + 3*x^2 + 4*x*y + 2*x + y^3 + 3*y^2 + 2*y

julia> char_poly = characteristic_polynomial(M)
q^3 - 6*q^2 + 11*q - 6

julia> factor(char_poly)
1 * (q - 2) * (q - 1) * (q - 3)

julia> A = chow_ring(M);

julia> GR, _ = graded_polynomial_ring(QQ,symbols(base_ring(A)));

julia> AA = chow_ring(M,ring=GR);

julia> vol_map = volume_map(M,AA)
#3605 (generic function with 1 method)

julia> e = matroid_groundset(M)[1];

julia> proper_flats = flats(M)[2:length(flats(M))-1];

julia> a = sum([AA[i] for i in 1:length(proper_flats) if e in proper_flats[i]])
x_{Edge(2, 1)} + x_{Edge(2, 1),Edge(3, 1),Edge(3, 2)} + x_{Edge(2, 1),Edge(4, 1),Edge(4, 2)} + x_{Edge(2, 1),Edge(4, 3)}

julia> b = sum([AA[i] for i in 1:length(proper_flats) if !(e in proper_flats[i])])
x_{Edge(3, 1)} + x_{Edge(3, 2)} + x_{Edge(4, 1)} + x_{Edge(4, 2)} + x_{Edge(4, 3)} + x_{Edge(3, 1),Edge(4, 2)} + x_{Edge(3, 1),Edge(4, 1),Edge(4, 3)} + x_{Edge(3, 2),Edge(4, 1)} + x_{Edge(3, 2),Edge(4, 2),Edge(4, 3)}

julia> k = 1
1

julia> R = base_ring(AA);

julia> g = grading_group(R)[1];

julia> PD1, mapPD1 = homogeneous_component(AA,k*g);

julia> basis_PD1 = [mapPD1(x) for x in gens(PD1)];

julia> PD2, mapPD2 = homogeneous_component(AA,(rank(M)-k-1)*g);

julia> basis_PD2 = [mapPD2(x) for x in gens(PD2)];

julia> Mat1 = matrix(QQ,[[vol_map(b1*b2) for b1 in basis_PD1] for b2 in basis_PD2])
[-1    0    0    0    0    0    0    1]
[ 0   -1    0    0    0    0    0    0]
[ 0    0   -1    0    0    0    0    1]
[ 0    0    0   -1    0    0    0    0]
[ 0    0    0    0   -1    0    0    1]
[ 0    0    0    0    0   -1    0    0]
[ 0    0    0    0    0    0   -1    0]
[ 1    0    1    0    1    0    0   -2]

julia> Mat2 = matrix(QQ,[[vol_map(b1*b^(rank(M)-2k-1)*b2) for b1 in basis_PD1] for b2 in basis_PD1]);

julia> Mat1 == Mat2
true

julia> RR, _ = graded_polynomial_ring(QQ,"y_#" => 1:length(basis_PD1));

julia> map = hom(RR,AA,basis_PD1);

julia> K = kernel(hom(RR,AA,[b^(rank(M)-2k)*b3 for b3 in basis_PD1]));

julia> basis_HR = [map(h) for h in gens(K) if degree(h).coeff==k*g.coeff]
7-element Vector{MPolyQuoRingElem{MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}}}:
 x_{Edge(4, 3)}
 -x_{Edge(2, 1),Edge(3, 1),Edge(3, 2)} + x_{Edge(2, 1),Edge(4, 1),Edge(4, 2)}
 -x_{Edge(2, 1),Edge(3, 1),Edge(3, 2)} + 2*x_{Edge(2, 1),Edge(4, 3)}
 -x_{Edge(2, 1),Edge(3, 1),Edge(3, 2)} + 2*x_{Edge(3, 1),Edge(4, 2)}
 -x_{Edge(2, 1),Edge(3, 1),Edge(3, 2)} + x_{Edge(3, 1),Edge(4, 1),Edge(4, 3)}
 -x_{Edge(2, 1),Edge(3, 1),Edge(3, 2)} + 2*x_{Edge(3, 2),Edge(4, 1)}
 -x_{Edge(2, 1),Edge(3, 1),Edge(3, 2)} + x_{Edge(3, 2),Edge(4, 2),Edge(4, 3)}

julia> Mat3 = matrix(QQ,[[(-1)^k*vol_map(b1*b^(rank(M)-2k-1)*b2) for b1 in basis_HR] for b2 in basis_HR])
[ 2   0   -2   0   -1   0   -1]
[ 0   2    1   1    1   1    1]
[-2   1    5   1    1   1    1]
[ 0   1    1   5    1   1    1]
[-1   1    1   1    2   1    1]
[ 0   1    1   1    1   5    1]
[-1   1    1   1    1   1    2]

julia> is_positive_definite(matrix(ZZ,[ZZ(i) for i in Mat3]))
true

julia> reduced_characteristic_polynomial(M)
q^2 - 5*q + 6

julia> [vol_map(a^(rank(M)-j-1)*b^j) for j in range(0,rank(M)-1)]
3-element Vector{QQFieldElem}:
 1
 5
 6
```
