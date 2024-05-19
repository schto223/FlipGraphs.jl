# DeltaComplexes

## Modelling a triangulation on a closed surface

A Δ-complex is a representation of a triangulation on a closed surface.
To define a triangulation on a closed surface, it does not suffice to take into account vertices and edges. We will also need to take into account the triangular faces between them. Therefore the triangulations are modeled using their dual graph.

Vertices are triangular faces which in turn consist of three points, and three edges. These points and edges are not necessarily distinct.
Edges in the dual Graph(i.e. the Δ-complex) connect two triangular faces if they in turn share an edge. In order to avoid confusion between the edges of the triangulation and the edges in the dual graph, I will, hence forwards call the latter *dual edge*

```@docs
    DeltaComplex
    get_num_edges
    get_num_points
    get_num_trifaces
    euler_characteristic
    genus
    adjacency_matrix
    diameter
    createDeltaComplex
    flip!
    is_flippable
```