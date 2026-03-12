---
layout: page
title: Introduction
---

## Authors: Claus Fieker and Michael Joswig

This page contains the code samples from the corresponding chapter in the OSCAR book. You can access the full chapter [here](https://link.springer.com/chapter/10.1007/978-3-031-62127-7_1).

Interested readers might also find the article [OSCAR: An introduction (M. Horn)](https://fachgruppe-computeralgebra.de/data/CA-Rundbrief/car72.pdf) entertaining.

```julia
julia> R, p = residue_ring(ZZ, 3);

julia> typeof(R)
zzModRing

julia> F, p = residue_field(ZZ, 3);

julia> typeof(F)
FqField

julia> typeof(2)
Int64

julia> 2^100
0

julia> ZZ(2)^100
1267650600228229401496703205376

julia> 1/2
0.5

julia> 1//2
1//2

julia> ccall(:sin, Float64, (Float64,), 1.0)
0.8414709848078965

julia> u = Int64[1,2,3,4,5,6];

julia> ccall(:memset, Cvoid, (Ptr{Int}, Int, UInt), u, 0, length(u)*sizeof(Int64))

julia> show(u)
[0, 0, 0, 0, 0, 0]

julia> prod(1-1/a^2 for a = 1:100 if is_prime(a))
0.6090337253995164

julia> QQx, x = polynomial_ring(QQ, "x");

julia> typeof(x)
QQPolyRingElem

julia> Ry, y = polynomial_ring(QQ, ["y"]);

julia> typeof(y)
Vector{QQMPolyRingElem} (alias for Array{QQMPolyRingElem, 1})

julia> QQx, x = QQ["x"];

julia> function next(a)
          return a+1
       end;

julia> next(2)
3

julia> next(x)
x + 1

julia> next('a')
'b': ASCII/Unicode U+0062 (category Ll: Letter, lowercase)

julia> next(sin)
ERROR: MethodError: no method matching +(::typeof(sin), ::Int64)
[...]

julia> function next(a::typeof(sin))
         return x->a(x)+1
       end;

julia> next(sin)
#55 (generic function with 1 method)

julia> ans(1)
1.8414709848078965

julia> sin(1)+1
1.8414709848078965

julia> next(cos)
ERROR: MethodError: no method matching +(::typeof(cos), ::Int64)
[...]

julia> supertype(typeof(sin))
Function

julia> supertype(typeof(cos))
Function

julia> function next(a::Function)
         return x->a(x)+1
       end;

julia> a = next(next(tan))
#57 (generic function with 1 method)

julia> a(1)
3.5574077246549023

julia> tan(1)+2
3.5574077246549023
```

```julia
julia> @time class_group(quadratic_field(7)[1]);
  0.659554 seconds (916.73 k allocations: 60.414 MiB, 4.44% gc time, 94.30% compilation time)

julia> @time class_group(quadratic_field(11)[1]);
  0.003981 seconds (82.12 k allocations: 4.215 MiB)
```

```julia
julia> R, mR = quo(ZZ, 2)
(Integers modulo 2, Map: ZZ -> ZZ/(2))

julia> S, mS = quo(ZZ, 3)
(Integers modulo 3, Map: ZZ -> ZZ/(3))

julia> typeof(R(1))
zzModRingElem

julia> typeof(S(1))
zzModRingElem

julia> S(1) + R(1)
ERROR: Operations on distinct residue rings not supported
Stacktrace:
 [1] error(s::String)
   @ Base ./error.jl:35
 [2] check_parent
   @ ~/.julia/packages/Nemo/F0iQ2/src/flint/nmod.jl:22 [inlined]
 [3] +(x::zzModRingElem, y::zzModRingElem)
   @ Nemo ~/.julia/packages/Nemo/F0iQ2/src/flint/nmod.jl:138
[...]

julia> one(Int)
1

julia> one(Float64)
1.0

julia> [S(1) R(1)]
1×2 Matrix{zzModRingElem}:
 1  1

julia> dot(ans, ans)
ERROR: Operations on distinct residue rings not supported
Stacktrace:
[...]

julia> matrix(S, 1, 2, [1, 2])
[1   2]

julia> dot(ans, ans)
2

julia> det([1 2; 3 4])
-2.0

julia> det(BigInt[1 2; 3 4])
-2

julia> inv(BigInt[1 2; 3 4])
2×2 Matrix{BigFloat}:
 -2.0   1.0
  1.5  -0.5

julia> inv(BigInt[1 2; 0 1])
2×2 Matrix{BigFloat}:
 1.0  -2.0
 0.0   1.0

julia> inv(ZZ[1 2; 3 4])
ERROR: Matrix not invertible
[...]

julia> inv(QQ[1 2; 3 4])
[  -2       1]
[3//2   -1//2]
```

```julia
julia> Qx,x = QQ[:x]
(Univariate polynomial ring in x over QQ, x)

julia> QQx, xx = QQ[:x]
(Univariate polynomial ring in x over QQ, x)

julia> Qx === QQx
true

julia> x == xx
true

julia> Qx, x = polynomial_ring(QQ, :x, cached = false)
(Univariate polynomial ring in x over QQ, x)

julia> Qx == QQx
false

julia> x + xx
ERROR: Incompatible polynomial rings in polynomial operation
[...]
```
