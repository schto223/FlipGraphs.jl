# FlipGraph

!!! warning "Warning"
    The flip graph code uses holeyDeltaComplex in order to differentiate between different triangulations.
    I've noticed an oversight, which makes renders the code to build flip graphs to produce flip graphs where certain vertices have been falsely glued together.

Contrary to flipgraphs of planar triangulations like that of a convex polygon, the flipgraph of a closed surface is generally infinitely large.
Therefore, it is impossible to construct the whole flipgraph. However, one can construct a local image of the flipgraph.

```@docs
    FlipGraph
    nv(::FlipGraph)
    ne(::FlipGraph)
    vertices(::FlipGraph)
    edges(::FlipGraph)
    has_vertex(::FlipGraph, ::Integer)
    has_edge(::FlipGraph,::Edge)
    has_edge(::FlipGraph,::Integer, ::Integer)
    has_edge(::FlipGraph,::HoleyDeltaComplex, ::HoleyDeltaComplex)
    neighbors(::FlipGraph, ::Integer)
    diameter(::FlipGraph)
```

## Comparing Triangulations

In order to decide whether two triangulations of a surface are in fact homeomorphic to each other and therefore equivalent, we need to use [`HoleyDeltaComplex`](@ref)'s.
One problem of deciding, if two triangulations are equivalent, is that the naming of the vertices, edges and points is completely arbitrary.
Unless we are interested in the modular flipgraph, we can assume that the points are identical in two different instances of `HoleyDeltaComplex`. 
However, in the act of flipping, `DualEdge`s and `TriFace`s "move around". It is therefore possible, to obtain two `HoleyDeltaComplex`s representing the same triangulation of a surface but with different labeled `TriFace`s and `DualEdge`s. 
Checking every possible permutation of `TriFace` and `DualEdge` labelings is not an option, as the number of possibilities would blow up immeadiately. 

What we do instead, is trying to find a canonical labeling. This Module uses a version of *McKay's canonical graph labeling algorithm*[^1] to try and determin a unique labeling based on the relationship to other vertices and edges. In general, it is not always possible to determin a unique labeling.However, with this method, we can reduce the number of labelings to a manageable number.

The mcKay Methods each return a permutation vector `p` which can be interpreted as Cauchy's one-line notation for permutations. 
For example `p = [3,5,1,2,6,4]` would correspond to the following permutation:

```math
\centering
σ = \begin{pmatrix}
1 & 2 & 3 & 4 & 5 & 6\\
3 & 5 & 1 & 2 & 6 & 4
\end{pmatrix} = (1 3)(2 5 6 4)
```

```@docs
    mcKay_points(::HoleyDeltaComplex; ::Bool)
    mcKay_vertices(::HoleyDeltaComplex; ::Bool)
    mcKay_edges(::HoleyDeltaComplex; ::Bool)

    is_isomorph(::HoleyDeltaComplex, ::HoleyDeltaComplex)
    is_isomorph_to(::HoleyDeltaComplex, ::HoleyDeltaComplex)
```

## Constructing the FlipGraph
Contrary to flipgraphs of planar triangulations like that of a convex polygon, the flipgraph of a closed surface is generally infinitely large.
Therefore, it is impossible to construct the whole flipgraph. However, one can construct a local image of the flipgraph.

```@docs
    flip_graph(::HoleyDeltaComplex, ::Integer; ::Bool)
```


[^1]: Hartke, S.G., & Radcliffe, A.J. (2008). McKay ’ s Canonical Graph Labeling Algorithm.