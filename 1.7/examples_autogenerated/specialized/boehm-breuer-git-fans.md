---
layout: page
title: Algorithms for GIT-Fans of Affine Torus Actions
---

## Authors: Janko Böhm and Thomas Breuer

This page contains the code samples from the corresponding chapter in the OSCAR book. You can access the full chapter [here](https://link.springer.com/chapter/10.1007/978-3-031-62127-7_8).


```julia
julia> using Oscar.GITFans

julia> Q = [
            1  1   0   0   0
            1  0   1   1   0
            1  0   1   0   1
            1  0   0   1   1
            0  1   0   0  -1
            0  1   0  -1   0
            0  1  -1   0   0
            0  0   1   0   0
            0  0   0   1   0
            0  0   0   0   1
            ];

julia> n = nrows(Q);

julia> Qt, T = polynomial_ring(QQ, :T => 1:n);

julia> D = free_abelian_group(ncols(Q));

julia> w = [D(Q[i, :]) for i = 1:n];

julia> R, T = grade(Qt, w);

julia> a = ideal([
               T[5]*T[10] - T[6]*T[9] + T[7]*T[8],
               T[1]*T[9]  - T[2]*T[7] + T[4]*T[5],
               T[1]*T[8]  - T[2]*T[6] + T[3]*T[5],
               T[1]*T[10] - T[3]*T[7] + T[4]*T[6],
               T[2]*T[10] - T[3]*T[9] + T[4]*T[8],
           ]);
```

```julia
julia> perms_list = [ [1,3,2,4,6,5,7,8,10,9], [5,7,1,6,9,2,8,4,10,3] ];

julia> sym = symmetric_group(n);

julia> G, emb = sub([sym(x) for x in perms_list]...);

julia> G
Permutation group of degree 10

julia> describe(G)
"S5"

julia> iso = isomorphism(G, symmetric_group(5));

julia> [iso(x) for x in gens(G)]
2-element Vector{PermGroupElem}:
 (1,3)
 (1,5,4,3,2)
```

```julia
julia> matrix_action = GITFans.action_on_target(Q, G)
Group homomorphism
  from permutation group of degree 10 and order 120
  to matrix group of degree 5 over QQ

julia> x = gen(G, 2)
(1,5,9,10,3)(2,7,8,4,6)

julia> Ax = matrix_action(x)
[-1   1   -1   -1   -2]
[ 1   0    1    1    1]
[ 1   0    0    1    1]
[ 0   0    0    0    1]
[ 1   0    1    0    1]

julia> matrix(QQ, Q[Vector{Int}(x), 1:5]) == matrix(QQ, Q) * Ax
true
```

```julia
julia> collector_cones = GITFans.orbit_cones(a, Q, G);

julia> length(collector_cones)
4

julia> orbit_list = GITFans.orbit_cone_orbits(collector_cones, matrix_action);

julia> println(map(length, orbit_list))
[10, 15, 10, 1]
```

```julia
julia> perm_actions = GITFans.action_on_orbit_cone_orbits(orbit_list, matrix_action);

julia> perm_actions[1]
Group homomorphism
  from permutation group of degree 10 and order 120
  to permutation group of degree 10
```

```julia
julia> q_cone = positive_hull(Q)
Polyhedral cone in ambient dimension 5

julia> (hash_list, edges) = GITFans.fan_traversal(orbit_list, q_cone, perm_actions);

julia> length(hash_list)
6

julia> println(edges)
Set([[4, 6], [2, 3], [3, 5], [1, 2], [3, 4]])
```

```julia
julia> fanobj = GITFans.hashes_to_polyhedral_fan(orbit_list, hash_list, matrix_action)
Polyhedral fan in ambient dimension 5

julia> println(f_vector(fanobj))
ZZRingElem[20, 110, 240, 225, 76]

julia> c = cones(fanobj, 5)[1]
Polyhedral cone in ambient dimension 5

julia> rays(c)
5-element SubObjectIterator{RayVector{QQFieldElem}}:
 [0, 1, 0, 0, -1]
 [1, 0, 1, 1, 0]
 [1, 1, 0, 0, 0]
 [1, 1, 0, 1, 0]
 [1, 1, 1, 0, 0]

julia> dim(fanobj)
5

julia> n_rays(fanobj)
20

julia> n_maximal_cones(fanobj)
76

julia> n_cones(fanobj)
672

julia> is_pointed(fanobj)
true

julia> is_complete(fanobj)
false
```

```julia
julia> fanobj = GITFans.git_fan(a, Q, G)
Polyhedral fan in ambient dimension 5
```

```julia
julia> expanded = GITFans.orbits_of_maximal_GIT_cones(orbit_list, hash_list, matrix_action);

julia> orbit_lengths = map(length, expanded); println(orbit_lengths)
[1, 10, 30, 20, 10, 5]

julia> sum(orbit_lengths)
76

julia> maxcones = vcat( expanded... );

julia> full_edges = GITFans.edges_intersection_graph(maxcones, size(Q, 2) - 1);

julia> length(full_edges)
180
```

```julia
julia> full_graph = Graph{Undirected}(length(maxcones));

julia> [add_edge!(full_graph, e...) for e in full_edges];

julia> visualize(full_graph)

julia> orbit_graph = Graph{Undirected}(length(hash_list));

julia> [add_edge!(orbit_graph, e...) for e in edges];

julia> visualize(orbit_graph)
```
