# FlipGraph

Contrary to flipgraphs of planar triangulations like that of a convex polygon, the flipgraph of a closed surface is generally infinitely large.
Therefore, it is impossible to construct the whole flipgraph. Not only are they infinitely in size, but they also "grow" very fast. It is therefore more interesting and achievable, to look at the modular flip graphs. These are graphs whose vertices are isotopy classes of triangulations.

```@docs
    FlipGraph
    FGVertex
    nv(::FlipGraph)
    ne(::FlipGraph)
    vertices(::FlipGraph)
    edges(::FlipGraph)
    has_vertex(::FlipGraph, ::Integer)
    has_edge(::FlipGraph,::Edge)
    has_edge(::FlipGraph,::Integer, ::Integer)
    neighbors(::FlipGraph, ::Integer)
    diameter(::FlipGraph)
```

## Comparing Triangulations


One problem of deciding, if two triangulations are equivalent, is that the naming of the vertices, edges and points is completely arbitrary.
In the act of flipping, `DualEdge`s and `TriFace`s "move around". It is therefore possible, to obtain two `HoleyDeltaComplex`s representing the same triangulation of a surface but with different labeled `TriFace`s and `DualEdge`s. 
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
    mcKay_points(::DeltaComplex; ::Bool)
    mcKay_vertices(::DeltaComplex; ::Bool)
    mcKay_edges(::DeltaComplex; ::Bool)

    is_isomorphic(::DeltaComplex, ::DeltaComplex)
    is_isomorphic(candidate::FGVertexCandidate, fgv::FGVertex)
```

## Constructing the FlipGraph

```@docs
    flipgraph_modular
```


[^1]: Hartke, S.G., & Radcliffe, A.J. (2008). McKay ’ s Canonical Graph Labeling Algorithm.