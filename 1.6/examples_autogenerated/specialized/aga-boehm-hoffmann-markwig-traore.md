---
layout: page
title: Algorithms for Gromov-Witten Invariants of Elliptic Curves
---

## Authors: Janko Böhm, Firoozeh Dastur, Alain Hoffmann, Hannah Markwig and Ali Traore

This page contains the code samples from the corresponding chapter in the OSCAR book. You can access the full chapter [here](https://link.springer.com/chapter/10.1007/978-3-031-62127-7_6).

Interested readers are referred to [this article](https://arxiv.org/abs/2311.11381) for more
details. For yet more background, see [this article](https://arxiv.org/abs/2309.03741).

```julia
using Pkg                            
Pkg.add(url = "https://github.com/singular-gpispace/GromovWitten", rev="d7af772"; io=devnull)
using GromovWitten
```

```julia
julia> using GromovWitten

julia> G =  feynman_graph([(1,2),(1,2),(1,2)]);

julia> G = feynman_graph([(1,2),(1,2),(2,4),(1,3),(3,4),(3,4)]);

julia> G = feynman_graph([(1,2),(1,2),(1,3),(2,4),(3,4),(3,5),(4,6),(5,6),(5,6)]);

julia> G = feynman_graph([(1,2),(1,3),(1,4),(2,4),(2,3),(3,4)]);
```

```julia
julia> using GromovWitten

julia> G = feynman_graph([(1, 3),(1, 2),(1, 2),(2, 4),(3, 4),(3, 4)])
FeynmanGraph([(1, 3), (1, 2), (1, 2), (2, 4), (3, 4), (3, 4)])

julia> F = feynman_integral(G);

julia> a = [0, 2, 1, 0, 0, 1];

julia> feynman_integral_branch_type(F, a)
256*q[2]^4*q[3]^2*q[6]^2

julia> feynman_integral_degree(F, 3)
288*q[1]^6 + 32*q[1]^4*q[2]^2 + 32*q[1]^4*q[3]^2 + 32*q[1]^4*q[5]^2 + 32*q[1]^4*q[6]^2 + 8*q[1]^2*q[2]^2*q[5]^2 + 8*q[1]^2*q[2]^2*q[6]^2 + 8*q[1]^2*q[3]^2*q[5]^2 + 8*q[1]^2*q[3]^2*q[6]^2 + 24*q[2]^6 + 152*q[2]^4*q[3]^2 + 8*q[2]^4*q[5]^2 + 8*q[2]^4*q[6]^2 + 152*q[2]^2*q[3]^4 + 32*q[2]^2*q[3]^2*q[5]^2 + 32*q[2]^2*q[3]^2*q[6]^2 + 32*q[2]^2*q[4]^4 + 8*q[2]^2*q[4]^2*q[5]^2 + 8*q[2]^2*q[4]^2*q[6]^2 + 8*q[2]^2*q[5]^4 + 32*q[2]^2*q[5]^2*q[6]^2 + 8*q[2]^2*q[6]^4 + 24*q[3]^6 + 8*q[3]^4*q[5]^2 + 8*q[3]^4*q[6]^2 + 32*q[3]^2*q[4]^4 + 8*q[3]^2*q[4]^2*q[5]^2 + 8*q[3]^2*q[4]^2*q[6]^2 + 8*q[3]^2*q[5]^4 + 32*q[3]^2*q[5]^2*q[6]^2 + 8*q[3]^2*q[6]^4 + 288*q[4]^6 + 32*q[4]^4*q[5]^2 + 32*q[4]^4*q[6]^2 + 24*q[5]^6 + 152*q[5]^4*q[6]^2 + 152*q[5]^2*q[6]^4 + 24*q[6]^6

julia> f = feynman_integral_degree_sum(F, 3)
288*q[1]^6 + 32*q[1]^4*q[2]^2 + 32*q[1]^4*q[3]^2 + 32*q[1]^4*q[5]^2 + 32*q[1]^4*q[6]^2 + 8*q[1]^4 + 8*q[1]^2*q[2]^2*q[5]^2 + 8*q[1]^2*q[2]^2*q[6]^2 + 8*q[1]^2*q[3]^2*q[5]^2 + 8*q[1]^2*q[3]^2*q[6]^2 + 24*q[2]^6 + 152*q[2]^4*q[3]^2 + 8*q[2]^4*q[5]^2 + 8*q[2]^4*q[6]^2 + 152*q[2]^2*q[3]^4 + 32*q[2]^2*q[3]^2*q[5]^2 + 32*q[2]^2*q[3]^2*q[6]^2 + 8*q[2]^2*q[3]^2 + 32*q[2]^2*q[4]^4 + 8*q[2]^2*q[4]^2*q[5]^2 + 8*q[2]^2*q[4]^2*q[6]^2 + 8*q[2]^2*q[5]^4 + 32*q[2]^2*q[5]^2*q[6]^2 + 8*q[2]^2*q[6]^4 + 24*q[3]^6 + 8*q[3]^4*q[5]^2 + 8*q[3]^4*q[6]^2 + 32*q[3]^2*q[4]^4 + 8*q[3]^2*q[4]^2*q[5]^2 + 8*q[3]^2*q[4]^2*q[6]^2 + 8*q[3]^2*q[5]^4 + 32*q[3]^2*q[5]^2*q[6]^2 + 8*q[3]^2*q[6]^4 + 288*q[4]^6 + 32*q[4]^4*q[5]^2 + 32*q[4]^4*q[6]^2 + 8*q[4]^4 + 24*q[5]^6 + 152*q[5]^4*q[6]^2 + 152*q[5]^2*q[6]^4 + 8*q[5]^2*q[6]^2 + 24*q[6]^6

julia> substitute(f)
1792*q[1]^6 + 32*q[1]^4
```
