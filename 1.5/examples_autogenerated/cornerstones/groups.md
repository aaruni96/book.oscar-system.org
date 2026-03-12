---
layout: page
title: Group Theory
---

## Authors: Max Horn and Claus Fieker

This page contains the code samples from the corresponding chapter in the OSCAR book. You can access the full chapter [here](https://link.springer.com/chapter/10.1007/978-3-031-62127-7_4).

A preprint of this chapter is available [here](https://arxiv.org/abs/2404.05871).

```julia
import Pkg
Pkg.add(name="GenericCharacterTables", version="0.4"; io=devnull)
using GenericCharacterTables
```
### 1 Introduction

```julia
julia> G_perm = @permutation_group(5, (2,5)(3,4), (1,3)(4,5))
Permutation group of degree 5

julia> order(G_perm)
10

julia> describe(G_perm)
"D10"

julia> r = G_perm[1] * G_perm[2]
(1,3,5,2,4)

julia> g = perm([1,2,3,5,4])  # 'tabular' notation
(4,5)

julia> g in G_perm
false

julia> G_perm([2,1,5,4,3])  # 'tabular' notation with explicit parent
(1,2)(3,5)

julia> cperm([1,2],[4,5])   # cycle decomposition
(1,2)(4,5)

julia> @perm (1,2)(4,5)     # cycle notation via macro
(1,2)(4,5)

julia> K = algebraic_closure(QQ)
Algebraic closure of rational field

julia> e = one(K)
{a1: 1.00000}

julia> s, c = sinpi(2*e/5), cospi(2*e/5)
({a4: 0.951057}, {a2: 0.309017})

julia> mat_rot = matrix([ c -s ; s c ]);

julia> mat_sigma1 = matrix(K, [ -1 0 ; 0 1 ]);

julia> G_mat = matrix_group(mat_rot, mat_sigma1)
Matrix group of degree 2
  over algebraic closure of rational field

julia> is_isomorphic(G_mat, G_perm)
true

julia> p = K.([0,1]);  # coordinates of vertex 1, expressed over K

julia> orb = orbit(G_mat, *, p)
G-set of
  matrix group of degree 2 over QQBar
  with seeds [[{a1: 0}, {a1: 1.00000}]]

julia> pts = collect(orb)
5-element Vector{Vector{QQBarFieldElem}}:
 [{a1: 0}, {a1: 1.00000}]
 [{a4: 0.951057}, {a2: 0.309017}]
 [{a4: 0.587785}, {a2: -0.809017}]
 [{a4: -0.951057}, {a2: 0.309017}]
 [{a4: -0.587785}, {a2: -0.809017}]

julia> visualize(convex_hull(pts))

julia> R2 = vector_space(K, 2); # the "euclidean" plane over K

julia> sigma_1 = hom(R2, R2, [-R2[1], R2[2]])
Module homomorphism
  from vector space of dimension 2 over QQBar
  to vector space of dimension 2 over QQBar

julia> rot = hom(R2, R2, mat_rot);

julia> G_generic = generic_group(closure([rot, sigma_1], *), *)[1]
Generic group of order 10 with multiplication table

julia> permutation_group(G_generic)
Permutation group of degree 10

julia> is_isomorphic(ans, G_perm)
true

julia> G_fp = fp_group(G_perm)
Finitely presented group of order 10

julia> gens(G_fp)
2-element Vector{FPGroupElem}:
 F1
 F2

julia> relators(G_fp)
3-element Vector{FPGroupElem}:
 F1^2
 F1^-1*F2*F1*F2^-4
 F2^5

julia> F = free_group(2)
Free group of rank 2

julia> G_cox, _ = quo(F, [F[1]^2, F[2]^2, (F[1]*F[2])^5])
(Finitely presented group, Hom: F -> G_cox)

julia> is_isomorphic(G_cox, G_perm)
true

julia> F = free_group(:a, :b, :c); a,b,c = gens(F);

julia> G, _ = quo(F, [a^2, b^2, c^2, (a*b)^3, (a*c)^2, (b*c)^3])
(Finitely presented group, Hom: F -> G)

julia> H, _ = quo(F, [a^2, b^2, c^2, (a*b)^3, (a*c)^3, (b*c)^3])
(Finitely presented group, Hom: F -> H)

julia> is_finite(G)
true

julia> describe(G)
"S4"

julia> describe(H)
"a finitely presented group"

julia> H1, _ = derived_subgroup(H)
(Sub-finitely presented group, Hom: H1 -> H)

julia> H2, _ = derived_subgroup(H1)
(Sub-finitely presented group, Hom: H2 -> H1)

julia> describe(H2)
"Z x Z"

julia> describe(quo(H, H2)[1])
"(C3 x C3) : C2"

julia> dihedral_group(10)
Pc group of order 10

julia> dihedral_group(PermGroup, 10)
Permutation group of degree 5
```
### 2 Group actions, 3 Databases

```julia
julia> gset(alternating_group(4))
G-set of
  alternating group of degree 4
  with seeds 1:4

julia> G = dihedral_group(6);

julia> U = sub(G, [g for g in gens(G) if order(g) == 2])[1]
Sub-pc group of order 2

julia> r = right_cosets(G, U)
Right cosets of
  sub-pc group of order 2 in
  pc group of order 6

julia> acting_group(r)
Pc group of order 6

julia> collect(r)
3-element Vector{GroupCoset{PcGroup, SubPcGroup, PcGroupElem}}:
 Right coset of U with representative <identity> of ...
 Right coset of U with representative f2
 Right coset of U with representative f2^2

julia> action_function(r)
* (generic function with 1736 methods)

julia> permutation(r, G[1])
(2,3)

julia> phi = right_coset_action(G,U)
Group homomorphism
  from pc group of order 6
  to permutation group of degree 3 and order 6

julia> phi(G[1]), phi(G[2])
((2,3), (1,2,3))

julia> function optimal_transitive_perm_rep(G)
         is_natural_symmetric_group(G) && return hom(G,G,gens(G))
         is_natural_alternating_group(G) && return hom(G,G,gens(G))
         cand = []  # pairs (U,h) with U ≤ G and h a map G -> Sym(G/U)
         for C in subgroup_classes(G)
           U = representative(C)
           h = right_coset_action(G, U)
           is_injective(h) && push!(cand, (U, h))
         end
         return argmax(a -> order(a[1]), cand)[2]
       end;

julia> U = dihedral_group(8)
Pc group of order 8

julia> optimal_transitive_perm_rep(U)
Group homomorphism
  from pc group of order 8
  to permutation group of degree 4 and order 8

julia> isomorphism(PermGroup, U)
Group homomorphism
  from pc group of order 8
  to permutation group of degree 4 and order 8

julia> permutation_group(U)
Permutation group of degree 4 and order 8

julia> for g in all_transitive_groups(degree => 3:8, !is_primitive)
         h = image(optimal_transitive_perm_rep(g))[1]
         if degree(h) < degree(g)
           id = transitive_group_identification(g)
           id_new = transitive_group_identification(h)
           println(id => id_new)
         end
       end
(6, 2) => (3, 2)
(6, 4) => (4, 4)
(6, 7) => (4, 5)
(6, 8) => (4, 5)
(8, 4) => (4, 3)
(8, 13) => (6, 6)
(8, 14) => (4, 5)
(8, 24) => (6, 11)
```

### 4 Character theory

#### SL(2,5)
The characters of SL(2,5) ... and the representation of chapter 8
```julia
julia> G = SL(2, 5)
SL(2,5)

julia> T = character_table(G)
Character table of SL(2,5)

  2  3                 1                 1  3                  1                  1  1  1  2
  3  1                 .                 .  1                  .                  .  1  1  .
  5  1                 1                 1  1                  1                  1  .  .  .
 
    1a               10a               10b 2a                 5a                 5b 3a 6a 4a
 
X_1  1                 1                 1  1                  1                  1  1  1  1
X_2  2 z_5^3 + z_5^2 + 1    -z_5^3 - z_5^2 -2 -z_5^3 - z_5^2 - 1      z_5^3 + z_5^2 -1  1  .
X_3  2    -z_5^3 - z_5^2 z_5^3 + z_5^2 + 1 -2      z_5^3 + z_5^2 -z_5^3 - z_5^2 - 1 -1  1  .
X_4  3    -z_5^3 - z_5^2 z_5^3 + z_5^2 + 1  3     -z_5^3 - z_5^2  z_5^3 + z_5^2 + 1  .  . -1
X_5  3 z_5^3 + z_5^2 + 1    -z_5^3 - z_5^2  3  z_5^3 + z_5^2 + 1     -z_5^3 - z_5^2  .  . -1
X_6  4                -1                -1  4                 -1                 -1  1  1  .
X_7  4                 1                 1 -4                 -1                 -1  1 -1  .
X_8  5                 .                 .  5                  .                  . -1 -1  1
X_9  6                -1                -1 -6                  1                  1  .  .  .

julia> R = gmodule(T[end])
G-module for G acting on vector space of dimension 6 over abelian closure of QQ

julia> S = gmodule(CyclotomicField, R)
G-module for G acting on vector space of dimension 6 over cyclotomic field of order 5

julia> schur_index(T[end])
2

julia> gmodule_minimal_field(S)
G-module for G acting on vector space of dimension 6 over number field

julia> B, mB = relative_brauer_group(base_ring(S), character_field(S));

julia> B
Relative Brauer group for cyclotomic field of order 5 over number field of degree 1 over QQ

julia> b = B(S)
Element of relative Brauer group of number field of degree 1 over QQ
  <2> -> 1//2 + Z
  <5> -> 0 + Z
  Real embedding of number field -> 1//2 + Z

julia> K = grunwald_wang(b)
Class field
  over number field with defining polynomial _$ - 1
    over rational field
with modulus
  finite part <4>
  infinite part
    [Infinite place of real embedding with 1.00]
with structure
  Z/2

julia> F, _ = absolute_simple_field(number_field(K))
(Number field of degree 2 over QQ, Map: F -> relative non-simple number field)

julia> L, _, F_to_L = compositum(base_ring(S), F)
(Number field of degree 8 over QQ, Map: cyclotomic field of order 5 -> L, Map: F -> L)

julia> gmodule_over(F_to_L, gmodule(L, S))
G-module for G acting on vector space of dimension 6 over F
```
#### Character theory

```julia
julia> t = character_table("J2");

julia> mx = character_table.(maxes(t));

julia> pis = [trivial_character(s)^t for s in mx];

julia> ords = orders_class_representatives(t);  println(ords)
[1, 2, 2, 3, 3, 4, 5, 5, 5, 5, 6, 6, 7, 8, 10, 10, 10, 10, 12, 15, 15]

julia> ord7 = findall(is_equal(7), ords);

julia> filt = filter(pi -> pi[ord7[1]] != 0, pis);

julia> degree.(filt)
2-element Vector{QQFieldElem}:
 100
 1800

julia> ord5 = findall(is_equal(5), ords);

julia> s = sum(filt);

julia> all(i -> s[i] == 0, ord5)
true

julia> (x -> [x["1a"], x["2b"], x["3a"], x["3b"]]).(filt)
2-element Vector{Vector{QQAbFieldElem{AbsSimpleNumFieldElem}}}:
 [100, 0, 10, 4]
 [1800, 20, 0, 6]

julia> ord3 = findall(is_equal(3), ords);

julia> [class_multiplication_coefficient(t, 3, i, ord7[1]) for i in ord3]
2-element Vector{ZZRingElem}:
 0
 70

julia> s = mx[findfirst(is_equal(filt[2]), pis)];

julia> identifier(s)
"L3(2).2"

julia> possible_class_fusions(s, t)
1-element Vector{Vector{Int64}}:
 [1, 2, 5, 6, 13, 3, 12, 14, 14]

julia> class_multiplication_coefficient(s, 6, 3, 5)
0
```

### 5 Applying finitely presented groups towards group extension

```julia
julia> Oscar.randseed!(42)

julia> G = dihedral_group(PermGroup, 10);

julia> h = epimorphism_from_free_group(G)
Group homomorphism
  from free group of rank 2
  to permutation group of degree 5

julia> gens(G)
2-element Vector{PermGroupElem}:
 (1,2,3,4,5)
 (2,5)(3,4)

julia> g = rand(G)
(1,4)(2,3)

julia> w = preimage(h, g)
x1^-1*x2^-1*x1^-3

julia> map_word(w, gens(G))  # evaluate w at generators to get back g
(1,4)(2,3)

julia> k, mk = kernel(h);  # mk is a map from k to domain(h)

julia> ngens(k)
11

julia> mk(k[2])            # pick any kernel generator ...
x1*x2^-2*x1^-1

julia> G[1]*G[-2]^2*G[-1]  # ... and it evaluates to the identity
()

julia> U = cyclic_group(5); V = cyclic_group(2); F = free_group(2);

julia> q = quo(F, [F[1]^5, F[2]^2*F[1]])[1];

julia> describe(q)
"C10"

julia> q, mq = quo(F, [F[1]^5, F[2]^2, F[-1]*F[-2]*F[1]*F[2]*F[1]^2])
(Finitely presented group, Hom: F -> q)

julia> describe(q)
"D10"
```
### 6 Computing cohomology

```julia
julia> irreducible_modules(symmetric_group(3))
3-element Vector{GModule}:
 G-module for symmetric group of degree 3 acting on vector space of dimension 1 over abelian closure of QQ
 G-module for symmetric group of degree 3 acting on vector space of dimension 1 over abelian closure of QQ
 G-module for symmetric group of degree 3 acting on vector space of dimension 2 over abelian closure of QQ

julia> G = symmetric_group(3);

julia> A = abelian_group([2,2,2])
(Z/2)^3

julia> M = gmodule(G, [hom(A, A, permuted(gens(A), g)) for g in gens(G)])
G-module for G acting on A
```

### 7 Finding irreducible representations of solvable groups

```julia
julia> function regular_gmodule(F::Field, G::Oscar.GAPGroup)
           FM = free_module(F, order(Int, G));
           gs = gset(G, *, [one(G)])
           return gmodule(G, [hom(FM, FM, permutation_matrix(F, permutation(gs, g))) for g in gens(G)])
       end;

julia> Base.:^(a::T, b::T) where T <: MatElem = inv(b)*a*b

julia> G = dihedral_group(10);

julia> M = regular_gmodule(GF(7), G);

julia> C = composition_factors_with_multiplicity(M)
3-element Vector{Any}:
 (G-module for G acting on vector space of dimension 1 over GF(7), 1)
 (G-module for G acting on vector space of dimension 1 over GF(7), 1)
 (G-module for G acting on vector space of dimension 4 over GF(7), 2)

julia> [is_absolutely_irreducible(x[1]) for x in C]
3-element Vector{Bool}:
 1
 1
 0

julia> phi = embed(GF(7), splitting_field(C[3][1]))
Morphism of finite fields
  from prime field of characteristic 7
  to finite field of degree 2 and characteristic 7

julia> M = extension_of_scalars(C[3][1], phi)
G-module for G acting on vector space of dimension 4 over GF(7, 2)

julia> composition_factors_with_multiplicity(M)
2-element Vector{Any}:
 (G-module for G acting on vector space of dimension 2 over GF(7, 2), 1)
 (G-module for G acting on vector space of dimension 2 over GF(7, 2), 1)

julia> G = pc_group(symmetric_group(4));

julia> C = abelian_closure(QQ)[1];

julia> F = free_module(C, 1);

julia> s, ms = sub(G, [G[3], G[4]]);   # subgroup generated by c, d

julia> T = trivial_gmodule(s, F);

julia> z = root_of_unity(C, 3);

julia> ss, mss = sub(G, [G[2], G[3], G[4]]);

julia> M = gmodule(ss, [hom(F, F, [z*F[1]]), action(T, s[1]), action(T, s[2])])
G-module for ss acting on F

julia> FF = free_module(C, 2);

julia> zm = 0*matrix(action(M, one(ss)));

julia> im = [hom(FF, FF, [matrix(action(M, x)) zm; zm matrix(action(M, preimage(mss, mss(x)^G[1])))]) for x in gens(ss)];

julia> pushfirst!(im, hom(FF, FF, matrix(C, [0 1; 1 0])));

julia> phi = gmodule(G, im);

julia> character(phi)
class_function(character table of G, [2, 0, -1, 2, 0])

julia> schur_index(ans)
1

julia> T = matrix(C, [C(1) -z; z+1 z]);

julia> [matrix(x)^T for x in action(phi)]
4-element Vector{AbstractAlgebra.Generic.MatSpaceElem{QQAbFieldElem{AbsSimpleNumFieldElem}}}:
 [1 0; -1 -1]
 [0 1; -1 -1]
 [1 0; 0 1]
 [1 0; 0 1]
```
### 8 Example: SL(2,5)

```julia
julia> G = SL(2, 5)
SL(2,5)

julia> T = character_table(G)
Character table of SL(2,5)

  2  3                 1                 1  3                  1                  1  1  1  2
  3  1                 .                 .  1                  .                  .  1  1  .
  5  1                 1                 1  1                  1                  1  .  .  .
 
    1a               10a               10b 2a                 5a                 5b 3a 6a 4a
 
X_1  1                 1                 1  1                  1                  1  1  1  1
X_2  2 z_5^3 + z_5^2 + 1    -z_5^3 - z_5^2 -2 -z_5^3 - z_5^2 - 1      z_5^3 + z_5^2 -1  1  .
X_3  2    -z_5^3 - z_5^2 z_5^3 + z_5^2 + 1 -2      z_5^3 + z_5^2 -z_5^3 - z_5^2 - 1 -1  1  .
X_4  3    -z_5^3 - z_5^2 z_5^3 + z_5^2 + 1  3     -z_5^3 - z_5^2  z_5^3 + z_5^2 + 1  .  . -1
X_5  3 z_5^3 + z_5^2 + 1    -z_5^3 - z_5^2  3  z_5^3 + z_5^2 + 1     -z_5^3 - z_5^2  .  . -1
X_6  4                -1                -1  4                 -1                 -1  1  1  .
X_7  4                 1                 1 -4                 -1                 -1  1 -1  .
X_8  5                 .                 .  5                  .                  . -1 -1  1
X_9  6                -1                -1 -6                  1                  1  .  .  .

julia> R = gmodule(T[end])
G-module for G acting on vector space of dimension 6 over abelian closure of QQ

julia> S = gmodule(CyclotomicField, R)
G-module for G acting on vector space of dimension 6 over cyclotomic field of order 5

julia> schur_index(T[end])
2

julia> gmodule_minimal_field(S)
G-module for G acting on vector space of dimension 6 over number field

julia> B, mB = relative_brauer_group(base_ring(S), character_field(S));

julia> B
Relative Brauer group for cyclotomic field of order 5 over number field of degree 1 over QQ

julia> b = B(S)
Element of relative Brauer group of number field of degree 1 over QQ
  <2> -> 1//2 + Z
  <5> -> 0 + Z
  Real embedding of number field -> 1//2 + Z

julia> K = grunwald_wang(b)
Class field
  over number field with defining polynomial _$ - 1
    over rational field
with modulus
  finite part <4>
  infinite part
    [Infinite place of real embedding with 1.00]
with structure
  Z/2

julia> F, _ = absolute_simple_field(number_field(K))
(Number field of degree 2 over QQ, Map: F -> relative non-simple number field)

julia> L, _, F_to_L = compositum(base_ring(S), F)
(Number field of degree 8 over QQ, Map: cyclotomic field of order 5 -> L, Map: F -> L)

julia> gmodule_over(F_to_L, gmodule(L, S))
G-module for G acting on vector space of dimension 6 over F
```

### 9 Generic character tables

This section showcases the
[`GenericCharacterTables.jl` package](https://github.com/oscar-system/GenericCharacterTables.jl).
An extended and possibly updated version of this example can be found in
[the package manual](https://oscar-system.github.io/GenericCharacterTables.jl/stable/book/).

```julia
julia> T = generic_character_table("SL3.n1")
Generic character table SL3.n1
  of order q^8 - q^6 - q^5 + q^3
  with 8 irreducible character types
  with 8 class types
  with parameters (a, b, m, n)

julia> T[4,4]
(q + 1)*exp(2π𝑖((a*n)//(q - 1))) + exp(2π𝑖((-2*a*n)//(q - 1)))

julia> h = T[2] * T[2];

julia> scalar_product(T[4], h)
0
With exceptions:
  2*n1 ∈ (q - 1)ℤ

julia> for i in 1:8 println("<$i, h> = ", scalar_product(T[i], h)) end
<1, h> = 1
<2, h> = 2
<3, h> = 2
<4, h> = 0
With exceptions:
  2*n1 ∈ (q - 1)ℤ
<5, h> = 0
With exceptions:
  2*n1 ∈ (q - 1)ℤ
<6, h> = 0
With exceptions:
  2*m1 - n1 ∈ (q - 1)ℤ
  m1 - 2*n1 ∈ (q - 1)ℤ
  m1 + n1 ∈ (q - 1)ℤ
  m1 ∈ (q - 1)ℤ
  m1 - n1 ∈ (q - 1)ℤ
  n1 ∈ (q - 1)ℤ
<7, h> = 0
With exceptions:
  n1 ∈ (q - 1)ℤ
<8, h> = 0
With exceptions:
  q*n1 ∈ (q^2 + q + 1)ℤ
  n1 ∈ (q^2 + q + 1)ℤ
  q*n1 + n1 ∈ (q^2 + q + 1)ℤ

julia> degree(linear_combination([1,2,2],[T[1],T[2],T[3]]))
2*q^3 + 2*q^2 + 2*q + 1

julia> degree(h)
q^4 + 2*q^3 + q^2

julia> parameters(T[4])
n ∈ {1,…, q - 1} except n ∈ (q - 1)ℤ

julia> T2 = set_congruence(T; remainder=0, modulus=2);

julia> (q, (a, b, m, n)) = parameters(T2);

julia> x = parameter(T2, "x");  # create an additional "free" variable

julia> s = specialize(T2[6], m, -n + (q-1)*x);  # force m = -n (mod q-1)

julia> scalar_product(s, T2(h))
1
With exceptions:
  3*n1 ∈ (q - 1)ℤ
  2*n1 ∈ (q - 1)ℤ
```
