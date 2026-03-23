---
layout: page
title: Polyhedral Geometry

---

## Authors: Taylor Brysiewicz and Michael Joswig

This page contains the code samples from the corresponding chapter in the OSCAR book. You can access the full chapter [here](https://link.springer.com/chapter/10.1007/978-3-031-62127-7_5).

A preprint of this chapter is available [here](https://arxiv.org/abs/2308.07459).

```julia
julia> points = [1 0; 1 1; 0 1; -1 0; -1 -1];

julia> P = convex_hull(points)
Polyhedron in ambient dimension 2

julia> facets(P)
5-element SubObjectIterator{AffineHalfspace{QQFieldElem}} over the halfspaces of R^2 described by:
-x_1 <= 1
-x_1 + x_2 <= 1
x_1 - 2*x_2 <= 1
x_1 <= 1
x_2 <= 1


julia> f_vector(P)
2-element Vector{ZZRingElem}:
 5
 5

julia> volume(P)
5//2

julia> lattice_points(P)
6-element SubObjectIterator{PointVector{ZZRingElem}}:
 [-1, -1]
 [-1, 0]
 [0, 0]
 [0, 1]
 [1, 0]
 [1, 1]
```

```julia
julia> P = polyhedron([1 5; 2 -1; -1 0; 0 -1], [20,10,0,0])
Polyhedron in ambient dimension 2

julia> LP = linear_program(P, [1,1], convention=:max)
Linear program
   max{c*x + k | x in P}
where P is a Polyhedron{QQFieldElem} and
   c=QQFieldElem[1, 1]
   k=0

julia> optimal_value(LP)
100//11

julia> optimal_vertex(LP)
2-element PointVector{QQFieldElem}:
 70//11
 30//11
```

```julia
julia> T = tetrahedron()
Polytope in ambient dimension 3

julia> vertices(T)
4-element SubObjectIterator{PointVector{QQFieldElem}}:
 [1, -1, -1]
 [-1, 1, -1]
 [-1, -1, 1]
 [1, 1, 1]
```

```julia
julia> [ f_vector(P) for P in [T, cube(3), cross_polytope(3),
                               dodecahedron(), icosahedron()] ]
5-element Vector{Vector{ZZRingElem}}:
 [4, 6, 4]
 [8, 12, 6]
 [6, 12, 8]
 [20, 30, 12]
 [12, 30, 20]
```

```julia
julia> D = dodecahedron()
Polytope in ambient dimension 3 with EmbeddedAbsSimpleNumFieldElem type coefficients

julia> vertices(D)[1]
3-element PointVector{EmbeddedAbsSimpleNumFieldElem}:
 1//2 (0.50)
 1//4*sqrt(5) + 3//4 (1.31)
 0 (0.00)
```

```julia
orbit_counts = MSet{Int}();

for k in 1:92
  J = johnson_solid(k)
  Aut = automorphism_group(J)[:on_vertices]
  push!(orbit_counts, length(orbits(Aut)))
end
```

```julia
julia> sort(collect(orbit_counts.dict))
14-element Vector{Pair{Int64, Int64}}:
  2 => 26
  3 => 24
  4 => 17
  5 => 8
  7 => 5
  8 => 2
  9 => 2
 12 => 1
 14 => 1
 15 => 1
 17 => 1
 20 => 1
 27 => 1
 29 => 2
```

```julia
julia> P3 = orbit_polytope([1,2,3,4], symmetric_group(4))
Polyhedron in ambient dimension 4

julia> dim(P3)
3

julia> describe(combinatorial_symmetries(P3))
"C2 x S4"
```

```julia
julia> Q = convex_hull([0 0 0], [0 1 0; 0 0 1], [1 0 0])
Polyhedron in ambient dimension 3

julia> vertices(Q)
0-element SubObjectIterator{PointVector{QQFieldElem}}

julia> minimal_faces(Q)
(base_points = PointVector{QQFieldElem}[[0, 0, 0]], lineality_basis = RayVector{QQFieldElem}[[1, 0, 0]])
```

```julia
julia> Rays = [1 0; 1 1; 0 1; -1 0; 0 -1];

julia> Cones = IncidenceMatrix([[1,2], [2,3], [3,4], [4,5], [5,1]]);

julia> PF = polyhedral_fan(Cones, Rays)
Polyhedral fan in ambient dimension 2
```

```julia
julia> P = rand_spherical_polytope(6,30)
Polytope in ambient dimension 6

julia> show(f_vector(P))
ZZRingElem[30, 336, 1468, 2874, 2568, 856]
julia> show(h_vector(P))
ZZRingElem[1, 24, 201, 404, 201, 24, 1]
julia> show(g_vector(P))
ZZRingElem[1, 23, 177, 203]
```

```julia
n_vertices = 30;
n_samples = 100;
g_vectors = Array{Int32}(undef,n_samples,2);

for i=1:n_samples
    RS = rand_spherical_polytope(6,n_vertices)
    g = g_vector(RS)
    g_vectors[i,1] = g[3] # notice index shift as Julia counts from 1
    g_vectors[i,2] = g[4]
end

using Plots
scatter(g_vectors[:,1], g_vectors[:,2],
        xlabel="g_2", ylabel="g_3", legend=false);

```

```julia
julia> min_g2 = minimum(g_vectors[:,1])
159

julia> max_g2 = maximum(g_vectors[:,1])
192

julia> show(upper_bound_g_vector(6,30))
[1, 23, 276, 2300]
julia> ub = [ Int(Polymake.polytope.pseudopower(g2,2)) for g2 in min_g2:max_g2 ]
34-element Vector{Int64}:
  990
  997
 1005
 1014
 1024
 1035
 1047
 1060
 1074
 1089
 1105
 1122
 1140
 1141
 1143
 1146
 1150
 1155
 1161
 1168
 1176
 1185
 1195
 1206
 1218
 1231
 1245
 1260
 1276
 1293
 1311
 1330
 1331
 1333
```

```julia
julia> lambda = Partition([3,1,1])
[3, 1, 1]

julia> GT = gelfand_tsetlin_polytope(lambda)
Polyhedron in ambient dimension 6

julia> lattice_points(GT)
6-element SubObjectIterator{PointVector{ZZRingElem}}:
 [3, 1, 1, 1, 1, 1]
 [3, 1, 1, 2, 1, 1]
 [3, 1, 1, 2, 1, 2]
 [3, 1, 1, 3, 1, 1]
 [3, 1, 1, 3, 1, 2]
 [3, 1, 1, 3, 1, 3]

julia> ehrhart_polynomial(GT)
2*x^2 + 3*x + 1

julia> volume(project_full(GT))
2
```

```julia
julia> w0 = @perm (1,3,2)
(1,3,2)

julia> genGT = gelfand_tsetlin_polytope(lambda,w0)
Polyhedron in ambient dimension 6

julia> lattice_points(genGT)
3-element SubObjectIterator{PointVector{ZZRingElem}}:
 [3, 1, 1, 3, 1, 1]
 [3, 1, 1, 3, 1, 2]
 [3, 1, 1, 3, 1, 3]
```

```julia
julia> dc = demazure_character(lambda,w0)
x1^3*x2*x3 + x1^2*x2^2*x3 + x1*x2^3*x3

julia> dc(1,1,1)
3
```

```julia
julia> C = cube(3, 0, 1);

julia> R, x, c = polynomial_ring(QQ, :x=>1:3, :c=>(0:1,0:1,0:1));

julia> f = sum(prod(x[i]^Int(p[i]) for i=1:3)
               * c[(Vector{Int}(p)+[1,1,1])...] for p=lattice_points(C))
x[1]*x[2]*x[3]*c[1, 1, 1] + x[1]*x[2]*c[1, 1, 0] + x[1]*x[3]*c[1, 0, 1] + x[1]*c[1, 0, 0] + x[2]*x[3]*c[0, 1, 1] + x[2]*c[0, 1, 0] + x[3]*c[0, 0, 1] + c[0, 0, 0]

julia> I = ideal(R, vcat([derivative(f,t) for t = x], [f]));

julia> D222 = eliminate(I,x)[1]
c[0, 0, 0]^2*c[1, 1, 1]^2 - 2*c[0, 0, 0]*c[1, 0, 0]*c[0, 1, 1]*c[1, 1, 1] - 2*c[0, 0, 0]*c[0, 1, 0]*c[1, 0, 1]*c[1, 1, 1] - 2*c[0, 0, 0]*c[1, 1, 0]*c[0, 0, 1]*c[1, 1, 1] + 4*c[0, 0, 0]*c[1, 1, 0]*c[1, 0, 1]*c[0, 1, 1] + c[1, 0, 0]^2*c[0, 1, 1]^2 + 4*c[1, 0, 0]*c[0, 1, 0]*c[0, 0, 1]*c[1, 1, 1] - 2*c[1, 0, 0]*c[0, 1, 0]*c[1, 0, 1]*c[0, 1, 1] - 2*c[1, 0, 0]*c[1, 1, 0]*c[0, 0, 1]*c[0, 1, 1] + c[0, 1, 0]^2*c[1, 0, 1]^2 - 2*c[0, 1, 0]*c[1, 1, 0]*c[0, 0, 1]*c[1, 0, 1] + c[1, 1, 0]^2*c[0, 0, 1]^2
```

```julia
julia> CubeFacets = [lattice_points(F) for F in faces(C,2)];

julia> VariableIndices = [[c+[1,1,1] for c=Vector{Vector{Int}}(F)] 
                           for F=CubeFacets];

julia> FacialDeterminants = [c[v[1]...]*c[v[4]...]-c[v[2]...]*c[v[3]...] 
                             for v=VariableIndices]
6-element Vector{QQMPolyRingElem}:
 c[0, 0, 0]*c[0, 1, 1] - c[0, 1, 0]*c[0, 0, 1]
 c[1, 0, 0]*c[1, 1, 1] - c[1, 1, 0]*c[1, 0, 1]
 c[0, 0, 0]*c[1, 0, 1] - c[1, 0, 0]*c[0, 0, 1]
 c[0, 1, 0]*c[1, 1, 1] - c[1, 1, 0]*c[0, 1, 1]
 c[0, 0, 0]*c[1, 1, 0] - c[1, 0, 0]*c[0, 1, 0]
 c[0, 0, 1]*c[1, 1, 1] - c[1, 0, 1]*c[0, 1, 1]

julia> E222 = D222*prod(FacialDeterminants);

julia> NP = newton_polytope(E222);

julia> dim(NP)
4

julia> show(f_vector(NP))
ZZRingElem[74, 152, 100, 22]
```

```julia
julia> C = cube(3,0,1)
Polytope in ambient dimension 3

julia> SP = secondary_polytope(C)
Polytope in ambient dimension 8

julia> V = point_matrix(vertices(NP))[:, 4:11].+1;

julia> SP == convex_hull(V)
true

julia> AllTriangulations = all_triangulations(C);

julia> Tri_as_SOP = [subdivision_of_points(C,T)
                     for T in AllTriangulations];

julia> GKZ_Vectors = [gkz_vector(T) for T in Tri_as_SOP];

julia> SP_from_GKZ = convex_hull(GKZ_Vectors)
Polyhedron in ambient dimension 8

julia> SP_from_GKZ == SP
true

julia> T = AllTriangulations[1]
6-element Vector{Vector{Int64}}:
 [1, 2, 3, 5]
 [2, 3, 4, 5]
 [2, 4, 5, 6]
 [3, 4, 5, 7]
 [4, 5, 6, 7]
 [4, 6, 7, 8]

julia> S = subdivision_of_points(C,T)
Subdivision of points in ambient dimension 3

julia> is_regular(S)
true

julia> show(min_weights(S))
[3, 1, 1, 0, 0, 0, 0, 1]
```

```julia
julia> IM = IncidenceMatrix(T);

julia> V = matrix(QQ, vertices(C));

julia> PC = polyhedral_complex(IM,V);

julia> visualize(PC)
```

```julia
julia> G = automorphism_group(C)[:on_vertices]
Permutation group of degree 8

julia> OrbitRepresentatives=
        unique([minimum(gset(G,permuted,[v])) 
          for v in GKZ_Vectors])
6-element Vector{Vector{QQFieldElem}}:
 [1, 3, 3, 5, 5, 3, 3, 1]
 [1, 3, 3, 5, 6, 2, 2, 2]
 [1, 3, 4, 4, 6, 2, 1, 3]
 [1, 4, 4, 3, 6, 1, 1, 4]
 [1, 5, 5, 1, 5, 1, 1, 5]
 [2, 2, 2, 6, 6, 2, 2, 2]

julia> OrbitSizes = 
        [length(
          filter(x->minimum(gset(G,permuted,[x]))==u,
          GKZ_Vectors)) for u in OrbitRepresentatives];

julia> show(OrbitSizes)
[12, 24, 24, 8, 2, 4]
```

```julia
julia> P = rand_spherical_polytope(6, 500; seed=11)
Polytope in ambient dimension 6

julia> Polymake.prefer("ppl") do
         @time n_facets(P)
       end
205.427732 seconds (7.13 G allocations: 321.966 GiB, 0.68% gc time)
58163

julia> P = rand_spherical_polytope(6, 500; seed=11);

julia> Polymake.prefer("beneath_beyond") do
         @time n_facets(P)
       end
 61.090983 seconds (335.63 M allocations: 5.309 GiB, 0.08% gc time)
58163

julia> P = rand_spherical_polytope(6, 500; seed=11);

julia> Polymake.prefer("libnormaliz") do
         @time n_facets(P)
       end
 11.167154 seconds (19.12 M allocations: 811.301 MiB, 0.07% gc time)
58163
```
