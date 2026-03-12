---
layout: page
title: Number Theory
---

## Authors: Claus Fieker and Tommy Hofmann

This page contains the code samples from the corresponding chapter in the OSCAR book. You can access the full chapter [here](https://link.springer.com/chapter/10.1007/978-3-031-62127-7_3).

A preprint of this chapter is available [here](https://arxiv.org/abs/2404.06858).

### 2 Number Fields and Rings of Integers

Basic examples and invariants.

```julia
julia> Qx, x = QQ["x"];

julia> K, a = number_field(x^2 - 235, "a")
(Number field of degree 2 over QQ, a)

julia> a^2 - 235 # confirm that a is a root of x^2 - 235
0

julia> b = 2 + 1//3 * a
1//3*a + 2

julia> coordinates(b)
2-element Vector{QQFieldElem}:
 2
 1//3

julia> K(x^3 + x + 2) # map the polynomial x^3 + x + 2 to K
236*a + 2

julia> degree(K)
2

julia> trace(a)
0

julia> norm(a)
-235

julia> OK = ring_of_integers(K);

julia> basis(OK)
2-element Vector{AbsSimpleNumFieldOrderElem}:
 1
 a

julia> discriminant(OK)
940

julia> prime_ideals_over(OK, 7)
2-element Vector{AbsSimpleNumFieldOrderIdeal}:
 <7, a + 2>
 <7, a + 5>

julia> factor(change_coefficient_ring(GF(7), x^2 - 235))
(x + 2) * (x + 5)

julia> factor(a * OK)
Dict{AbsSimpleNumFieldOrderIdeal, Int64} with 2 entries:
  <47, a> => 1
  <5, a>  => 1

julia> A, m = class_group(OK);

julia> A
Z/6

julia> m(zero(A)) # apply m to the neutral element of A
Ideal of maximal order of number field of degree 2 over QQ
  of norm 1
  of minimum 1
with generator 1

julia> P = prime_ideals_over(OK, 2)[1]
Ideal of maximal order of number field of degree 2 over QQ
  of norm 2
  of minimum 2
with 2-normal generators [2, a + 1]

julia> preimage(m, P)
Abelian group element [3]

julia> is_principal_with_data(P^2)
(true, 2)

julia> P^2 == 2*OK
true

julia> U, mU = unit_group(OK);

julia> U
Z/2 x Z

julia> mU(U[1]), mU(U[2])
(-1, 3*a + 46)

julia> preimage(mU, -214841715*a - 3293461126)
Abelian group element [1, 5]

julia> -214841715*a - 3293461126 == (-1)^1 * (3a + 46)^5
true
```
#### Figure 1

Code to plot the lattice of integers.

```julia
using Plots;
w1, w2 = Complex(1),(sqrt(Complex(-3)) + 1)/2
pts = filter(z -> max(abs(imag(z)), abs(real(z))) <= 3,
             [ a*w1 + b*w2 for a in -5:5 for b in -5:5]);
x = range(0, real(w1) + real(w2), 100);
y1 = (x -> min(sqrt(3)*x, imag(w2))).(x);
y2 = (x -> max(0, sqrt(3)*x - sqrt(3))).(x);
p = scatter(real.(pts), imag.(pts), framestyle=:origin,
            legend = false, mc = :red)
plot!(p, x, y1, fillrange = y2, fillalpha = 0.25, color = :blue)
```
#### Figure 2, 3

Code to plot the unit group in a quadratic field: first the conjugates

```julia
using Plots;
pts = [ (a, b) for a in -8:8 for b in -5:5]
x = range(-7.5, 7.5, 500)
y1 = (x -> x^2 - 1 < 0 ? missing : sqrt((x^2 - 1)/3)).(x)
y2 = (x -> sqrt((x^2 + 1)/3)).(x)
scatter(pts, framestyle=:origin, leg = false, mc = :blue, alpha = 0.25)
plot!(x, [y1 -y1 y2 -y2], color = :black)
units = [(1, 0), (-1, 0), (2, 1), (-2, 1), (-2, -1),
         (2, -1), (7, 4), (7, -4), (-7, 4), (-7, -4)];
scatter!(units, color = :red)
```

Code to plot the unit group in a quadratic field: then the logarithm.

```julia
using Plots
l(a, b) = (log(abs(a + sqrt(3)*b)), log(abs(a - sqrt(3)*b)))
pts = filter!(z -> maximum(abs.(z)) < 3,
              [l(a, b) for a in -50:50 for b in -50:50])
scatter(pts, framestyle=:origin, leg = false, mc = :blue, alpha = 0.25)
x = range(-3, 3, 100); y = -x; plot!(x, y, color = :black)
units = [ i .* l(2, 1) for i in -2:2]
scatter!(units, color = :red)
```
### 3 General Number Fields and Morphisms

```julia
julia> Qx, x = QQ["x"];

julia> K, a = number_field([x^2 - 2, x^2 - 3]);

julia> a[1]^2 == 2 && a[2]^2 == 3
true

julia> k, b = quadratic_field(3);

julia> kt, t = k["t"];

julia> L, c = number_field(t^2 - 2);

julia> base_field(L) == k
true

julia> base_field(K) == QQ
true

julia> degree(K)
4

julia> degree(L)
2

julia> absolute_degree(K)
4

julia> absolute_norm(c)
4

julia> f = hom(K, K, [-a[1], -a[2]]);

julia> f(a[1]) == -a[1] && f(a[2]) == -a[2]
true

julia> g = hom(L, L, hom(k, k, -b), -c);

julia> g(c) == -c && g(L(b)) == L(-b)
true

julia> Lprime, LprimetoL = absolute_simple_field(L);

julia> defining_polynomial(Lprime)
x^4 - 10*x^2 + 1

julia> C, _ = class_group(ring_of_integers(Lprime)); order(C)
1
```

### 4 Embeddings of Number Fields

#### 4.1

```julia
julia> Qx, x = QQ["x"];

julia> K, a = number_field(x^3 - 2, "a");

julia> signature(K)
(1, 1)

julia> real_embeddings(K)
1-element Vector{AbsSimpleNumFieldEmbedding}:
 Real embedding with 1.26 of K

julia> embs = complex_embeddings(K) # printed with limited precision
3-element Vector{AbsSimpleNumFieldEmbedding}:
 Real embedding with 1.26 of K
 Imaginary embedding with -0.63 + 1.09 * i of K
 Imaginary embedding with -0.63 - 1.09 * i of K

julia> embs[1](a - 1)
[0.2599210499 +/- 8.44e-11]
```
#### 4.2

```julia
julia> Qx, x = QQ["x"];

julia> K, a = embedded_number_field(x^3 - 2, 1.26);

julia> 0 < a < 2
true
```

```julia
julia> Qx, x = QQ["x"];

julia> K, a = number_field(x^4 - 13*x^2 + 16, "a");

julia> G, m = automorphism_group(PermGroup, K); G # we only print G
Permutation group of degree 4

julia> describe(G)
"C2 x C2"

julia> phi = m(G[1])
Map
  from number field of degree 4 over QQ
  to number field of degree 4 over QQ
defined by
  a -> -a

julia> phi(a)
-a

julia> preimage(m, hom(K, K, -a))
(1,3)(2,4)

julia> k, a = number_field(x^2 - 18, "a"); kt, t = k["t"];

julia> K, b = number_field(t^4 + (a + 6)*t^2 + 2a + 9, "b");

julia> G, m = automorphism_group(PermGroup, K); describe(G)
"C4"

julia> G, m = absolute_automorphism_group(PermGroup, K); describe(G)
"Q8"

julia> K, a = number_field(x^4 - 2);

julia> G, mA = automorphism_group(PermGroup, K);

julia> describe(G)
"C2"

julia> G, C = galois_group(K); describe(G)
"D8"

julia> L = splitting_field(x^4 - 2);

julia> GL, = automorphism_group(PermGroup, L);

julia> fl, = is_isomorphic(G, GL); fl
true

julia> N, = normal_closure(K);

julia> GN, = automorphism_group(PermGroup, N);

julia> fl, = is_isomorphic(G, GN); fl
true

julia> rts = roots(C, 2)
4-element Vector{QadicFieldElem}:
 (8*11^0 + O(11^2))*a + 8*11^0 + 9*11^1 + O(11^2)
 (3*11^0 + 10*11^1 + O(11^2))*a + 7*11^0 + 4*11^1 + O(11^2)
 (3*11^0 + 10*11^1 + O(11^2))*a + 3*11^0 + 11^1 + O(11^2)
 (8*11^0 + O(11^2))*a + 4*11^0 + 6*11^1 + O(11^2)

julia> parent(rts[1])
Unramified extension of 11-adic numbers of degree 2

julia> k, = number_field(x^8 + x^7 - 7*x^6 - 23*x^5 + 22*x^4 + 80*x^3 + 99*x^2 + 27*x + 8, "a"); # not appearing in book

julia> defining_polynomial(k)
x^8 + x^7 - 7*x^6 - 23*x^5 + 22*x^4 + 80*x^3 + 99*x^2 + 27*x + 8

julia> describe(galois_group(k)[1])
"SL(2,3)"
```
### 5 Symmetries

```julia
Qx, x = QQ["x"]
K, a = number_field(x^9 - 3*x^8 + x^6 + 15*x^5 - 13*x^4 -
                    3*x^3 + 4*x - 1, "a")
G, C = galois_group(K)
subsG = subgroups(G)
H = first(H for H in subsG if order(H) == 27)
k, = simplify(fixed_field(C, H))
```

### 6 Examples
In the following it is necessary to have the package `Tally` installed, e.g. by running `]add Tally`.

#### 6.1
```julia
julia> fields = [d for d in 1:10^6 if
                 d != 1 && is_fundamental_discriminant(d)];

julia> length(fields)
303957

julia> clgrps = [class_group(quadratic_field(d, cached = false)[1])[1]
                 for d in fields];

julia> using Tally;

julia> tally(clgrps,
             by = x -> sylow_subgroup(x, 5)[1],
             equivalence = (x, y) -> is_isomorphic(x, y)[1])
Tally with 303957 items in 4 groups:
[Z/1]     | 291400 | 95.869%
[Z/5]     |  12324 |  4.055%
[Z/25]    |    230 |  0.076%
[(Z/5)^2] |      3 |  0.001%

julia> ab_grps = abelian_group.([[1], [5], [25], [5, 5]]);

julia> w = prod(1 - 1/5.0^i for i in 2:1000)
0.9504159948390403

julia> [ 100 * w/(1.0 * (order(A) * 
         order(automorphism_group(A)))) for A in ab_grps]
4-element Vector{BigFloat}:
 95.041[...]
  4.752[...]
  0.190[...]
  0.007[...]

julia> clnumbs = order.(clgrps);

julia> x = range(1, length(clnumbs), 1000);

julia> values = [count(c -> is_divisible_by(c, 5),
                       clnumbs[1:Int(ceil(step))])/Int(ceil(step)) for
                 step in x];

julia> pr = 1 - prod(1 - 1/5.0^i for i in 2:1000)
0.04958400516095973

julia> using Plots;
       plot(x, [values fill(pr, length(x))], linewidth = 2,
            xlabel = "Number of fields",
            label = ["Proportion" "Prediction"])
Plot{Plots.GRBackend() n=2}
```
#### 6.2

```julia
julia> Oscar.randseed!(VERSION > v"1.13.0-" ? 1145739 : 3371100);

julia> Qx, x = QQ["x"];

julia> K, a = number_field(x^4 + 4*x^2 + 2, "a");

julia> b = rand(K, -10:10) # random element with coefficients between -10 and 10
5*a^3 - 7*a^2 + 8*a - 8

julia> A, mA = automorphism_group(K);

julia> list = [mA(s)(b) for s in A]
4-element Vector{AbsSimpleNumFieldElem}:
 5*a^3 - 7*a^2 + 8*a - 8
 -7*a^3 + 7*a^2 - 26*a + 20
 -5*a^3 - 7*a^2 - 8*a - 8
 7*a^3 + 7*a^2 + 26*a + 20

julia> X = matrix(QQ, coordinates.(list))
[-8     8   -7    5]
[20   -26    7   -7]
[-8    -8   -7   -5]
[20    26    7    7]

julia> det(X) != 0
true

julia> det(matrix(QQ, [coordinates(mA(s)(a)) for s in A]))
0

julia> V, f = galois_module(K);

julia> OK = ring_of_integers(K);

julia> M = f(OK);

julia> is_free(M)
false

julia> fl = is_locally_free(M, 2)
false

julia> prime_decomposition_type(OK, 2)
1-element Vector{Tuple{Int64, Int64}}:
 (1, 4)

julia> K, a = number_field(x^4 - x^3 + x^2 - x + 1, "a");

julia> is_tamely_ramified(K)
true

julia> V, f = galois_module(K); OK = ring_of_integers(K); M = f(OK);

julia> fl, c = is_free_with_basis(M); # the elements of c form a basis

julia> b = preimage(f, c[1]) # the element may be different per session
a^3

julia> A, mA = automorphism_group(K);

julia> X = matrix(ZZ, [coordinates(OK(mA(s)(b))) for s in A])
[0    0    0    1]
[1   -1    1   -1]
[0    0   -1    0]
[0    1    0    0]

julia> det(X)
-1

julia> k, = number_field(x^4 - 13*x^2 + 16); candidates = []; for F in abelian_normal_extensions(k, [2], ZZ(10)^13)
         K, = absolute_simple_field(number_field(F))
         if !is_tamely_ramified(K)
           continue
         end
         if !is_isomorphic(galois_group(K)[1], quaternion_group(8))
           continue
         end
         push!(candidates, K)
       end

julia> length(candidates)
2

julia> K = candidates[2]; V, f = galois_module(K); OK = ring_of_integers(K); M = f(OK);

julia> fl = is_free(M)
false

julia> defining_polynomial(K)
x^8 + 105*x^6 + 3465*x^4 + 44100*x^2 + 176400
```

```julia
for i in 1:8, j in 1:number_of_small_groups(i)
  G = small_group(i, j)
  if is_abelian(G)
    continue
  end
  QG = QQ[G]
  ZG = integral_group_ring(QG)
  println(i, " ", j, " ", describe(G), ": ", locally_free_class_group(ZG))
end

```

```julia
Qx, x = QQ["x"]
k, = number_field(x^4 - 13*x^2 + 16)
candidates = []
for F in abelian_normal_extensions(k, [2], ZZ(10)^13)
  K, = absolute_simple_field(number_field(F))
  if !is_tamely_ramified(K)
    continue
  end
  if !is_isomorphic(galois_group(K)[1], quaternion_group(8))
    continue
  end
  push!(candidates, K)
end
```
