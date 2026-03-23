---
layout: page
title: Tropical Implicitization Revisited

---

## Authors: Kemal Rose, Bernd Sturmfels and Simon Telen

This page contains the code samples from the corresponding chapter in the OSCAR book. You can access the full chapter [here](https://link.springer.com/chapter/10.1007/978-3-031-62127-7_17).

```julia
using Pkg
Pkg.add(url = "https://github.com/kemalrose/TropicalImplicitization.jl"; io=devnull)
using TropicalImplicitization
```

```julia
julia> using TropicalImplicitization, Oscar;

julia> R, (t,) = polynomial_ring(QQ,["t"]);

julia> f1 = 11*t^2 + 5*t^3 - 1*t^4;

julia> f2 = 11 + 11*t + 7*t^8;

julia> Q1 = newton_polytope(f1);

julia> Q2 = newton_polytope(f2);

julia> newton_pols = [Q1, Q2];

julia> cone_list, weight_list = get_tropical_cycle(newton_pols);

julia> Delta = get_polytope_from_cycle(cone_list, weight_list);

julia> B = lattice_points(Delta);

julia> n_samples = length(B)-1;

julia> P = sample([f1,f2], n_samples);

julia> M_BP = get_vandermonde_matrix(B,P);

julia> coeffs_F = nullspace(M_BP)[2];
```

```julia
julia> A = [1 1 1 1 1 1 1 1; 0 0 0 0 1 1 1 1; 
            0 0 1 1 0 0 1 1; 0 1 0 1 0 1 0 1];

julia> cone_list, weight_list = get_trop_A_disc(A);
```

```julia
julia> Oscar.randseed!(43);

julia> Delta = get_polytope_from_cycle(cone_list, weight_list);

julia> f_vec, lattice_pts = f_vector(Delta), lattice_points(Delta);

julia> print(f_vec);
ZZRingElem[6, 14, 16, 8]
```

```julia
julia> Delta = get_polytope_from_cycle(cone_list, weight_list);
```

```julia
julia> A = [1 1 1 1 1 1; 2 3 5 7 11 13; 13 8 5 3 2 1];

julia> cone_list, weight_list = get_trop_A_disc(A);

julia> Delta = get_polytope_from_cycle(cone_list, weight_list);

julia> @time mons, coeffs = compute_A_discriminant(A,Delta,GF(101));
Collect all lattice points.
Sample 2754.0 points from Discriminant
Construct Vandermonde matrix.
Compute coefficients of the Discriminant.
 97.247222 seconds (210.17 M allocations: 17.505 GiB, 6.51% gc time, 0.00% compilation time)
```

```julia
julia> verts1 = [898 -614; -570 817; 892 -594];

julia> verts2 = [-603 -481; -623 -127; -36 732];

julia> verts3 = [-548 -864; -151 873; 800 -861];

julia> pols = convex_hull.([verts1, verts2, verts3]);

julia> Delta = get_polytope_from_cycle(get_tropical_cycle(pols)...);

julia> print(f_vector(Delta))
ZZRingElem[25, 49, 26]
```

```julia
julia> cone_list = positive_hull.([[1, 1, 0],[1, 2, 3], [1,0,1], [-1, -1, -4//3]]);

julia> weight_list = ones(Int64, 4);

julia> cone_list, weight_list = get_chow_fan(cone_list,weight_list);
```

```julia
julia> C_translated = get_polytope_from_cycle(cone_list,weight_list);
```
