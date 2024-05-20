# HoleyDeltaComplex

The [`DeltaComplex`](@ref) structure is usefull for many applications. However, it has a slight issue that two different triangulations of the same surface, may result in the same DeltaComplex. In other words, a DeltaComplex misses some information about the underlying triangulation and can therefore not be used if one would, for example, like to construct the flipgraph (or rather a portion of latter) of the triangulations on a given surface.

The [`HoleyDeltaComplex`](@ref) structure solves this issue, by keeping track of the holes through which every single edge passes, and the orders thereof. This comes of course with a decrease in efficiency. It is therefore advised to use the simpler [`DeltaComplex`](@ref) if this bonus information is not needed. Another drawback is, that the `HoleyDeltaComplex` structure currently does not work for non-orientable surfaces.
As the `HoleyDeltaComplex` is essentially an extension of the `DeltaComplex` structure, all methods which take DeltaComplexes also work with HoleyDeltaComplexes  

```@docs
    HoleyDeltaComplex
    Hole
    Crossing
    holey_delta_complex

    subdivide!(::HoleyDeltaComplex, ::Integer)
    rename_edges!(::HoleyDeltaComplex)
    rename_points!(::HoleyDeltaComplex)
    rename_vertices!(::HoleyDeltaComplex)
    np(::HoleyDeltaComplex)
    nv(::HoleyDeltaComplex)
    ne(::HoleyDeltaComplex)
    points(::HoleyDeltaComplex)
    get_point(::HoleyDeltaComplex)
    vertices(::HoleyDeltaComplex)
    get_vertex(::HoleyDeltaComplex)
    edges(::HoleyDeltaComplex)
    get_edge(::HoleyDeltaComplex)
    edges_id(::HoleyDeltaComplex)
    get_edge_id(::HoleyDeltaComplex)
    id(::HoleyDeltaComplex)
    other_endpoint(::HoleyDeltaComplex)
    _is_similar(::HoleyDeltaComplex)

    euler_characteristic(::HoleyDeltaComplex)
    genus(::HoleyDeltaComplex)
    diameter_triangulation(::HoleyDeltaComplex)
    diameter_deltaComplex(::HoleyDeltaComplex)
    adjacency_matrix_deltaComplex(::HoleyDeltaComplex)
    adjacency_matrix_triangulation(::HoleyDeltaComplex)
    multi_adjacency_matrix_triangulation(::HoleyDeltaComplex)
    point_degrees(::HoleyDeltaComplex)
    relative_point_degrees(::HoleyDeltaComplex)

    flip!(::HoleyDeltaComplex)
    is_flippable(::HoleyDeltaComplex)
    random_flips!(::HoleyDeltaComplex)
    randomize!(::HoleyDeltaComplex)

    num_crossings(::HoleyDeltaComplex)
    edge_crossings(::HoleyDeltaComplex)
    get_crossing(::HoleyDeltaComplex)
    remove_holeloops!(::HoleyDeltaComplex)

    is_isomorph(::HoleyDeltaComplex)
    is_isomorph_to(::HoleyDeltaComplex)
```