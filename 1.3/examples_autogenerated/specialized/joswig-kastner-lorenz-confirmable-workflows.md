---
layout: page
title: Confirmable Workflows
---

## Authors: Michael Joswig, Lars Kastner and Benjamin Lorenz

This page contains the code samples from the corresponding chapter in the OSCAR book. You can access the full chapter [here](https://link.springer.com/chapter/10.1007/978-3-031-62127-7_19).

```julia
julia> Oscar.versioninfo(full=true)
OSCAR version 1.0.0
  combining:
    AbstractAlgebra.jl   v0.40.1
    GAP.jl               v0.10.3
    Hecke.jl             v0.30.2
    Nemo.jl              v0.43.1
    Polymake.jl          v0.11.14
    Singular.jl          v0.22.4
  building on:
    Antic_jll               v0.201.500+0
    Arb_jll                 v200.2300.0+0
    Calcium_jll             v0.401.100+0
    FLINT_jll               v200.900.9+0
    GAP_jll                 v400.1200.200+9
    Singular_jll            v403.214.1400+0
    libpolymake_julia_jll   v0.11.4+0
    libsingular_julia_jll   v0.43.0+0
    polymake_jll            v400.1100.1+0
See `]st -m` for a full list of dependencies.

Julia Version 1.10.2
Commit bd47eca2c8a (2024-03-01 10:14 UTC)
Build Info:
  Official https://julialang.org/ release
Platform Info:
  OS: Linux (x86_64-linux-gnu)
  CPU: 8 × 11th Gen Intel(R) Core(TM) i7-1165G7 @ 2.80GHz
  WORD_SIZE: 64
  LIBM: libopenlibm
  LLVM: libLLVM-15.0.7 (ORCJIT, tigerlake)
Threads: 1 default, 0 interactive, 1 GC (on 8 virtual cores)
Official https://julialang.org/ release
```

```julia
julia> F, o = finite_field(7,2)
(Finite field of degree 2 and characteristic 7, o)

julia> R, (y,z) = F["y","z"]
(Multivariate polynomial ring in 2 variables over F, FqMPolyRingElem[y, z])

julia> save("p.mrdi", 2*y^3*z^4 + (o + 3)*z^2 + 5*o*y + 1)

julia> save("q.mrdi", y^17*z^3 + (o + 4)*z + 2*o*y)

julia> print( load("p.mrdi") + load("q.mrdi") )
y^17*z^3 + 2*y^3*z^4 + (o + 3)*z^2 + (o + 4)*z + 1
```

```julia
julia> m = matrix(ZZ, [1 2; 3 4])
[1   2]
[3   4]

julia> snf(m)
[1   0]
[0   2]

julia> methods(snf, Tuple{MatrixElem{ZZRingElem}})
# 2 methods for generic function "snf" from AbstractAlgebra:
 [1] snf(x::ZZMatrix)
     @ Nemo ~/.julia/packages/Nemo/xwSoX/src/flint/fmpz_mat.jl:1096
 [2] snf(a::MatrixElem{T}) where T<:RingElement
     @ ~/.julia/packages/AbstractAlgebra/dsta0/src/Matrix.jl:5431

julia> invoke(snf, Tuple{MatrixElem{ZZRingElem}}, m) == snf(m)
true
```
