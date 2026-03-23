---
layout: page
title: Monomial Bases in Lie Theory

---

## Authors: Xin Fang and Ghislain Fourier

This page contains the code samples from the corresponding chapter in the OSCAR book. You can access the full chapter [here](https://link.springer.com/chapter/10.1007/978-3-031-62127-7_13).

The interested reader can find more information in [this article](https://arxiv.org/pdf/2403.15018).

```julia
julia> basis_lie_highest_weight_operators(:A, 4)
10-element Vector{Tuple{Int64, Vector{QQFieldElem}}}:
 (1, [1, 0, 0, 0])
 (2, [0, 1, 0, 0])
 (3, [0, 0, 1, 0])
 (4, [0, 0, 0, 1])
 (5, [1, 1, 0, 0])
 (6, [0, 1, 1, 0])
 (7, [0, 0, 1, 1])
 (8, [1, 1, 1, 0])
 (9, [0, 1, 1, 1])
 (10, [1, 1, 1, 1])
```

```julia
julia> basis_lie_highest_weight(:A, 4, [2,1,2,1], [1,2,3,4,1,5,8,2,6,3]; monomial_ordering=:degrevlex)
Monomial basis of a highest weight module
  of highest weight 2*w_1 + w_2 + 2*w_3 + w_4
  of dimension 8750
  with monomial ordering degrevlex([x1, x2, x3, x4, x5, x6, x7, x8, x9, x10])
over abstract Lie algebra of type A4 over QQ
where the used birational sequence consists of the following roots:
  [a_1, a_2, a_3, a_4, a_1, a_1 + a_2, a_1 + a_2 + a_3, a_2, a_2 + a_3, a_3]
  
```

```julia
julia> basis_lie_highest_weight_ffl(:A, 3, [1,1,1])
Monomial basis of a highest weight module
  of highest weight w_1 + w_2 + w_3
  of dimension 64
  with monomial ordering degrevlex([x1, x2, x3, x4, x5, x6])
over abstract Lie algebra of type A3 over QQ
where the used birational sequence consists of the following roots:
  [a_1 + a_2 + a_3, a_2 + a_3, a_1 + a_2, a_3, a_2, a_1]
```

```julia
julia> basis_lie_highest_weight_string(:B, 3, [1,1,1], [3,2,3,2,1,2,3,2,1])
Monomial basis of a highest weight module
  of highest weight w_1 + w_2 + w_3
  of dimension 512
  with monomial ordering neglex([x1, x2, x3, x4, x5, x6, x7, x8, x9])
over abstract Lie algebra of type B3 over QQ
where the used birational sequence consists of the following roots:
  [a_3, a_2, a_3, a_2, a_1, a_2, a_3, a_2, a_1]
```

```julia
julia> basis_lie_highest_weight_lusztig(:D, 4, [1,1,1,1], [4,3,2,4,3,2,1,2,4,3,2,1])
Monomial basis of a highest weight module
  of highest weight w_1 + w_2 + w_3 + w_4
  of dimension 4096
  with monomial ordering wdegrevlex([x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12], [1, 1, 3, 2, 2, 1, 5, 4, 3, 3, 2, 1])
over abstract Lie algebra of type D4 over QQ
where the used birational sequence consists of the following roots:
  [a_4, a_3, a_2 + a_3 + a_4, a_2 + a_3, a_2 + a_4, a_2, a_1 + 2*a_2 + a_3 + a_4, a_1 + a_2 + a_3 + a_4, a_1 + a_2 + a_4, a_1 + a_
  2 + a_3, a_1 + a_2, a_1]
```

```julia
julia> basis_lie_highest_weight_nz(:C, 3, [1,1,1], [3,2,3,2,1,2,3,2,1])
Monomial basis of a highest weight module
  of highest weight w_1 + w_2 + w_3
  of dimension 512
  with monomial ordering degrevlex([x1, x2, x3, x4, x5, x6, x7, x8, x9])
over abstract Lie algebra of type C3 over QQ
where the used birational sequence consists of the following roots:
  [a_3, a_2, a_3, a_2, a_1, a_2, a_3, a_2, a_1]
```
