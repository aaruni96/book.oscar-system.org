---
layout: page
title: Generic Root Counts of Tropically Transverse Systems
---

## Authors: Isaac Holt and Yue Ren

```julia
using Pkg
Pkg.add("MixedSubdivisions"; io=devnull)

import LibGit2
LibGit2.clone("https://github.com/isaacholt100/generic_root_count", "generic_root_count")
include(joinpath(@__DIR__,"generic_root_count","src","main.jl"));
```

```julia
julia> T = tropical_semiring()
Min tropical semiring

julia> T(1)+T(2), T(1)*T(2), zero(T), one(T)
((1), (3), infty, (0))

julia> T = tropical_semiring(max)
Max tropical semiring

julia> T(1)+T(2), T(1)*T(2), zero(T), one(T)
((2), (3), -infty, (0))
```

```julia
julia> T = tropical_semiring()
Min tropical semiring

julia> M = matrix(T,[0 1; 2 3])
[(0)   (1)]
[(2)   (3)]

julia> v = T.([-1,-1])
2-element Vector{TropicalSemiringElem{typeof(min)}}:
 (-1)
 (-1)

julia> M^2*v
2-element Vector{TropicalSemiringElem{typeof(min)}}:
 (-1)
 (1)

julia> R,(x,y) = T["x","y"];

julia> f = 1*x^2+2*y^2+0
(1)*x^2 + (2)*y^2 + (0)

julia> f^2+2*f
(2)*x^4 + (3)*x^2*y^2 + (4)*y^4 + (1)*x^2 + (2)*y^2 + (0)

julia> evaluate(f,T.([1,1]))
(0)
```

```julia
julia> T = tropical_semiring()
Min tropical semiring

julia> M = matrix(T,[0 1; 2 3])
[(0)   (1)]
[(2)   (3)]

julia> det(M)
(3)
```

```julia
julia> nu = tropical_semiring_map(GF(3))
Map into Min tropical semiring encoding the trivial valuation on Prime field of characteristic 3

julia> nu.([0,1,2])
3-element Vector{TropicalSemiringElem{typeof(min)}}:
 infty
 (0)
 (0)
```

```julia
julia> nu = tropical_semiring_map(QQ,3,max)
Map into Max tropical semiring encoding the 3-adic valuation on Rational field

julia> nu.([1//3,1,3])
3-element Vector{TropicalSemiringElem{typeof(max)}}:
 (1)
 (0)
 (-1)
```

```julia
julia> Ft,t = rational_function_field(GF(3),"t")
(Rational function field over GF(3), t)

julia> nu = tropical_semiring_map(Ft,t,max)
Map into Max tropical semiring encoding the t-adic valuation on Rational function field over GF(3)

julia> nu.([t^(-1),1,t])
3-element Vector{TropicalSemiringElem{typeof(max)}}:
 (1)
 (0)
 (-1)
```

```julia
julia> R,(x1,x2,x3,x4) = QQ["x1","x2","x3","x4"];

julia> I = ideal([x1+2*x2-3*x3, 3*x2-4*x3+5*x4]);

julia> nu_2 = tropical_semiring_map(QQ,2);

julia> GBI = groebner_basis(I,nu_2,[0,0,0,0])
2-element Vector{QQMPolyRingElem}:
 x1 + 2*x2 - 3*x3
 3*x2 - 4*x3 + 5*x4

julia> inI = initial(I,nu_2,[0,0,0,0])
Ideal generated by
  x1 + x3
  x2 + x4

julia> GBI = groebner_basis(I,nu_2,[1,0,0,1])
2-element Vector{QQMPolyRingElem}:
 3*x2 - 4*x3 + 5*x4
 -x1 - 2*x2 + 3*x3

julia> inI = initial(I,nu_2,[1,0,0,1])
Ideal generated by
  x2
  x3

julia> coefficient_ring(inI)
Prime field of characteristic 2
```

```julia
julia> K,t = rational_function_field(GF(101),"t");

julia> nu = tropical_semiring_map(K,t);

julia> R,(x,y,z) = K["x","y","z"];

julia> I = intersect(ideal([x+y+z+1,2*x+11*y+23*z+31]),
                     ideal([t^3*x*y*z-1]));

julia> TropV = tropical_variety(I,nu)
2-element Vector{TropicalVariety}:
 Min tropical variety
 Min tropical variety

julia> dim.(TropV)
2-element Vector{Int64}:
 2
 1
```

```julia
julia> K,t = rational_function_field(QQ,"t");

julia> nu = tropical_semiring_map(K,t,max);

julia> R,(x,y) = K["x","y"];

julia> f = t^3+x+t^2*y+x*(x^2+y^2)
x^3 + x*y^2 + x + t^2*y + t^3

julia> TropH = tropical_hypersurface(f,nu)
Max tropical hypersurface

julia> vertices_and_rays(TropH) # 1,2,3 are vertices, rest are rays
7-element SubObjectIterator{Union{PointVector{QQFieldElem}, RayVector{QQFieldElem}}}:
 [0, -1]
 [-2, 0]
 [-3, -1]
 [0, 0]
 [1, 1]
 [-1, 1]
 [-1, 0]

julia> maximal_polyhedra(IncidenceMatrix,TropH)
7×7 IncidenceMatrix
[4, 5]
[1, 4]
[2, 4]
[2, 6]
[2, 3]
[1, 3]
[3, 7]
```

```julia
julia> tropf = tropical_polynomial(f,nu)
x^3 + x*y^2 + x + (-2)*y + (-3)

julia> tropical_hypersurface(tropf)
Max tropical hypersurface
```

```julia
julia> points = matrix(ZZ,collect(exponents(tropf)))
[3   0]
[1   2]
[1   0]
[0   1]
[0   0]

julia> heights = QQ.(collect(coefficients(tropf)))
5-element Vector{QQFieldElem}:
 0
 0
 0
 -2
 -3

julia> Delta = subdivision_of_points(points,heights)
Subdivision of points in ambient dimension 2

julia> tropical_hypersurface(Delta,max)
Max tropical hypersurface
```

```julia
julia> R,(w,x,y,z) = QQ["w","x","y","z"];

julia> nu = tropical_semiring_map(QQ);

julia> I = ideal([w+x+y+z,w+2*x+5*y+11*z]);

julia> TropL = tropical_linear_space(I,nu)
Min tropical linear space

julia> vertices_and_rays(TropL)
5-element SubObjectIterator{Union{PointVector{QQFieldElem}, RayVector{QQFieldElem}}}:
 [0, 0, 0, 0]
 [0, -1, -1, -1]
 [0, 1, 0, 0]
 [0, 0, 1, 0]
 [0, 0, 0, 1]

julia> maximal_polyhedra(IncidenceMatrix,TropL)
4×5 IncidenceMatrix
[1, 2]
[1, 3]
[1, 4]
[1, 5]
```

```julia
julia> A = matrix(QQ,[1 1 1 0; 1 1 0 1])
[1   1   1   0]
[1   1   0   1]

julia> TropL1 = tropical_linear_space(A, nu) # not Stiefel
Min tropical linear space

julia> tropA = nu.(A)
[(0)   (0)     (0)   infty]
[(0)   (0)   infty     (0)]

julia> TropL2 = tropical_linear_space(tropA) # not equal TropL1
Min tropical linear space
```

```julia
julia> plueckerIndices = [[1,2],[1,3],[1,4],[2,3],[2,4],[3,4]];

julia> algebraicPlueckerVector = [0,-1,-1,1,1,1];

julia> TropL1 = tropical_linear_space(plueckerIndices,
                                      algebraicPlueckerVector,
                                      nu)
Min tropical linear space

julia> tropicalPlueckerVector = nu.(algebraicPlueckerVector);

julia> TropL2 = tropical_linear_space(plueckerIndices,
                                      tropicalPlueckerVector)
Min tropical linear space
```

```julia
julia> R,(x,y) = QQ["x","y"];

julia> nu = tropical_semiring_map(QQ);

julia> f = 1+x+y+x*(x^2+y^2);

julia> TropH = tropical_hypersurface(f,nu)
Min tropical hypersurface

julia> TropC = tropical_curve(TropH)
Embedded min tropical curve
```

```julia
julia> G = dualgraph(cube(3))
Undirected graph with 6 nodes and the following edges:
(3, 1)(3, 2)(4, 1)(4, 2)(5, 1)(5, 2)(5, 3)(5, 4)(6, 1)(6, 2)(6, 3)(6, 4)

julia> tropical_curve(G)
Abstract min tropical curve
```

```julia
julia> K,t = rational_function_field(QQ,"t");

julia> nu = tropical_semiring_map(K,t,max);

julia> R,(x,y) = K["x","y"];

julia> f = t^3+x+t^2*y+x*(x^2+y^2)
x^3 + x*y^2 + x + t^2*y + t^3

julia> g = t^4+t^4*x+t^2*y+y*(x^2+y^2)
x^2*y + t^4*x + y^3 + t^2*y + t^4

julia> TropHf = tropical_hypersurface(f,nu)
Max tropical hypersurface

julia> TropHg = tropical_hypersurface(g,nu)
Max tropical hypersurface

julia> TropV = stable_intersection(TropHf,TropHg)
Max tropical variety

julia> vertices(TropV)
4-element SubObjectIterator{PointVector{QQFieldElem}}:
 [0, 0]
 [0, -4]
 [-3, -1]
 [-3, -2]

julia> multiplicities(TropV) # same ordering as above
4-element Vector{ZZRingElem}:
 4
 2
 2
 1
```

```julia
julia> include("generic_root_count/src/main.jl")

julia> A,(a0,a1,a2,a3,a4,b0,b1,b2,b3,b4) = QQ["a0","a1","a2","a3","a4","b0","b1","b2","b3","b4"];

julia> R,(u,v) = A["u","v"];

julia> f = a0+a1*u+a2*v+a3*u*(u^2+v^2)+a4*u*(u^2+v^2)^2;

julia> g = b0+b1*u+b2*v+b3*v*(u^2+v^2)+b4*v*(u^2+v^2)^2;

julia> generic_root_count([f,g])
9
```
