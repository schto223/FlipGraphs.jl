# FlipGraph

Contrary to flipgraphs of planar triangulations like that of a convex polygon, the flipgraph of a closed surface is generally infinitely large.
Therefore, it is impossible to construct the whole flipgraph. However, one can construct a local image of the flipgraph.

```@docs
    FlipGraph
    flip_graph
    vertices
    edges
    has_edge
    diameter
    mcKay_points
    mcKay_vertices
    mcKay_edges
    invert_perm
```