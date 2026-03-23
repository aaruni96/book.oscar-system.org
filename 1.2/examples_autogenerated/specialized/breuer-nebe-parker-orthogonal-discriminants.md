---
layout: page
title: An Atlas of Orthogonal Representations

---

## Authors: Thomas Breuer, Gabriele Nebe and Richard Parker

This page contains the code samples from the corresponding chapter in the OSCAR book. You can access the full chapter [here](https://link.springer.com/chapter/10.1007/978-3-031-62127-7_10).

```julia
julia> ch = character_table("Co2", 2)[2];

julia> degree(ch)
22

julia> Oscar.OrthogonalDiscriminants.od_from_order(ch)
(true, "O+")
```

```julia
julia> ch = character_table("Co3", 3)[2];

julia> degree(ch)
22

julia> Oscar.OrthogonalDiscriminants.od_from_eigenvalues(ch)
(true, "O+")
```

```julia
julia> ch = character_table("A12")[26];

julia> degree(ch)
1728

julia> Oscar.OrthogonalDiscriminants.od_for_specht_module(ch)
(true, "1")
```

```julia
julia> ch = character_table("R(27)")[16];

julia> degree(ch)
18278

julia> Oscar.OrthogonalDiscriminants.od_from_p_subgroup(ch, 3)
(true, "-3")
```

```julia
julia> Oscar.OrthogonalDiscriminants.show_with_ODs(
       character_table("G2(3)", 2))
G2(3)mod2

            2   6   3   3   .  1  1  .  .  .  .   .   .
            3   6   6   6   6  4  4  .  3  3  3   .   .
            7   1   .   .   .  .  .  1  .  .  .   .   .
           13   1   .   .   .  .  .  .  .  .  .   1   1
             
               1a  3a  3b  3c 3d 3e 7a 9a 9b 9c 13a 13b
           2P  1a  3a  3b  3c 3d 3e 7a 9a 9c 9b 13b 13a
           3P  1a  1a  1a  1a 1a 1a 7a 3c 3c 3c 13a 13b
           7P  1a  3a  3b  3c 3d 3e 1a 9a 9b 9c 13b 13a
          13P  1a  3a  3b  3c 3d 3e 7a 9a 9b 9c  1a  1a
     d OD   2                                          
 X_1 1      +   1   1   1   1  1  1  1  1  1  1   1   1
 X_2 1 O-   +  14   5   5  -4  2 -1  .  2 -1 -1   1   1
 X_3 2      o  64  -8  -8   1  4 -2  1  1  A /A  -1  -1
 X_4 2      o  64  -8  -8   1  4 -2  1  1 /A  A  -1  -1
 X_5 1 O-   +  78  -3  -3  -3 -3  6  1  .  .  .   .   .
 X_6 1 O+   +  90   9   9   9  .  . -1  .  .  .  -1  -1
 X_7 1 O-   +  90  -9  18   .  3 -3 -1 -3  .  .  -1  -1
 X_8 1 O-   +  90  18  -9   .  3 -3 -1 -3  .  .  -1  -1
 X_9 1 O-   + 378  -9  -9   9 -3 -6  .  3  .  .   1   1
X_10 2 O+   + 448  16  16 -11 -2 -2  .  1  1  1   B  B*
X_11 2 O+   + 448  16  16 -11 -2 -2  .  1  1  1  B*   B
X_12 1 O+   + 832 -32 -32  -5  4  4 -1  1  1  1   .   .

A = 3z_3 + 1
/A = -3z_3 - 2
B = -z_13^11 - z_13^8 - z_13^7 - z_13^6 - z_13^5 - z_13^2 - 1
B* = z_13^11 + z_13^8 + z_13^7 + z_13^6 + z_13^5 + z_13^2
```

```julia
julia> show_OD_info("G2(3)")
G2(3):  2^6*3^6*7*13
--------------------

 i| chi|     K|disc|           2|3|            7|      13
--+----+------+----+------------+-+-------------+--------
 2| 14a|     Q|  -3|         14a| |          14a|     14a
  |    |      |    |          O-| |           O+|      O+
--+----+------+----+------------+-+-------------+--------
 5| 78a|     Q|  -3|         78a| |          78a|     78a
  |    |      |    |          O-| |           O+|      O+
--+----+------+----+------------+-+-------------+--------
 9|104a|     Q|  21|     14a+90a| |     (def. 1)|    104a
  |    |      |    |      O-, O+| |             |      O-
--+----+------+----+------------+-+-------------+--------
10|168a|     Q|  13|     78a+90a| |         168a|(def. 1)
  |    |      |    |      O-, O+| |           O-|        
--+----+------+----+------------+-+-------------+--------
11|182a|     Q|  -3| 14a+78a+90c| |         182a|    182a
  |    |      |    |  O-, O-, O-| |           O+|      O+
--+----+------+----+------------+-+-------------+--------
12|182b|     Q|  -3| 14a+78a+90b| |         182b|    182b
  |    |      |    |  O-, O-, O-| |           O+|      O+
--+----+------+----+------------+-+-------------+--------
15|448a|Q(b13)|   1|        448a| |         448a|14a+434a
  |    |      |    |          O+| |           O+|  O+, O+
--+----+------+----+------------+-+-------------+--------
16|448b|Q(b13)|   1|        448b| |         448b|14a+434a
  |    |      |    |          O+| |           O+|  O+, O+
--+----+------+----+------------+-+-------------+--------
17|546a|     Q|  -3|78a+90b+378a| |         546a|    546a
  |    |      |    |  O-, O-, O-| |           O+|      O+
--+----+------+----+------------+-+-------------+--------
18|546b|     Q|  -3|78a+90c+378a| |         546b|    546b
  |    |      |    |  O-, O-, O-| |           O+|      O+
--+----+------+----+------------+-+-------------+--------
19|728a|     Q|   1|    14a+378a| |         728a|    728a
  |    |      |    |      O-, O-| |           O+|      O+
--+----+------+----+------------+-+-------------+--------
20|728b|     Q|   1|    14a+378a| |         728b|    728b
  |    |      |    |      O-, O-| |           O+|      O+
--+----+------+----+------------+-+-------------+--------
23|832a|     Q|   1|        832a| |64ab+78a+626a|    832a
  |    |      |    |          O+| |   O+, O+, O+|      O+
```

```julia
julia> K, _ = quadratic_field(13)
(Real quadratic field defined by x^2 - 13, sqrt(13))

julia> ray_class_field(3*maximal_order(K))
Class field defined mod (<3, 3>, InfPlc{AbsSimpleNumField, AbsSimpleNumFieldEmbedding}[]) of structure Z/1
```

```julia
julia> function is_galois_discriminant_field(data)
         chi = Oscar.OrthogonalDiscriminants.character_of_entry(data)
         F, emb = character_field(chi)
         c = conductor(emb(gen(F)))
         galgens = Oscar.AbelianClosure.generators_galois_group_cyclotomic_field(c)
         delta = atlas_irrationality(data[:valuestring])
         return all(x -> is_square(preimage(emb, delta * x(delta))),
                    galgens)
       end;

julia> info = all_od_infos(characteristic => 0, is_simple);

julia> filter!(r -> r[:valuestring] != "?" &&
                    conductor(atlas_irrationality(r[:valuestring])) > 1,
               info);

julia> length(info)
58

julia> filter!(!is_galois_discriminant_field, info);

julia> length(info)
26

julia> println(sort!(collect(Set([r[:groupname] for r in info]))))
["HN", "He", "J1", "J3", "ON", "Ru"]
```

```julia
julia> info = all_od_infos(characteristic => 0, degree => 1);

julia> all(x -> x[:valuestring] == "?" ||
                is_odd(parse(Int, x[:valuestring])),
           info)
true
```

```julia
julia> plus = [];  minus = [];

julia> for d in all_od_infos()
         if d[:valuestring] == "O+"
           push!(plus, (d[:groupname], d[:characteristic], d[:degree],
                        parse(Int, filter(isdigit, d[:charname]))))
         elseif d[:valuestring] == "O-"
           push!(minus, (d[:groupname], d[:characteristic], d[:degree],
                         parse(Int, filter(isdigit, d[:charname]))))
         end
       end

julia> both = intersect!(plus, minus);

julia> filter(x -> x[2] == 2, both)
1-element Vector{Any}:
 ("G2(3)", 2, 1, 90)

julia> length(both)
103
```
