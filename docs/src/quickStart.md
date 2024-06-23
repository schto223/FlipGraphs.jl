# Quick Start

```@meta
DocTestSetup = quote
    using FlipGraphs
end
```

If you're already familiar with the concept of **flip graphs**, **triangulations** on **closed surfaces**, and **Δ-complexes**, and don't want to read the whole documentation, then here are some quick examples of what you can do with this package.

In any other case, please be sure to have a look at the rest of the documentation first.

## Triangulated Convex Polygon

Create a triangulated convex 10-gon:

```jldoctest ggg
julia> g = triangulated_polygon(8)
TriangulatedPolygon with 8 vertices, and adjacency list:
 1  → [2, 8]
 2  → [1, 3, 8]
 3  → [2, 4, 8, 7]
 4  → [3, 5, 7, 6]
 5  → [4, 6]
 6  → [5, 7, 4]
 7  → [6, 8, 3, 4]
 8  → [7, 1, 2, 3]
```

Check if the edge going from vertex 2 to vertex 10 can be flipped:

```jldoctest ggg
julia> is_flippable(g, 2, 8)
true
```

Flip the edge connecting vertices 2 and 10:

```jldoctest ggg
julia> flip!(g, 2, 8)
TriangulatedPolygon with 8 vertices, and adjacency list:
 1  → [2, 8, 3]
 2  → [1, 3]
 3  → [2, 4, 8, 7, 1]
 4  → [3, 5, 7, 6]
 5  → [4, 6]
 6  → [5, 7, 4]
 7  → [6, 8, 3, 4]
 8  → [7, 1, 3]
```

Construct the flip graph of a convex octagon:

```jldoctest
julia> G = flipgraph_planar(8)
FlipGraphPlanar with 132 vertices and 330 edges
```

Export the generated flip graph as a .gml file:

```julia-repl
julia> export_gml("C:/Users/USERNAME/Desktop/FILENAME.gml", G);
```


## Δ-Complex / Triangulation of closed surface

### `DeltaComplex`

A `DeltaComplex` is the dual of a triangulation on a closed surface.
It can be used to compute things like the diameter, but it does not offer a unique model of a triangulation on a closed surface. 
Every DeltaComplex can be interpreted as the homeomorphism class of triangulations of points on a closed surface.

Create a `DeltaComplex` of a surface of genus 1 with 2 points:

```jldoctest DDD
julia> D = deltacomplex(1, 2)
DeltaComplex on orientable surface of genus 1 with 2 points
4 TriFaces:
 TriFace #1: Points(1 1 2) Neighbors(2 3 4)
 TriFace #2: Points(1 1 1) Neighbors(4 1 3)
 TriFace #3: Points(1 1 2) Neighbors(2 4 1)
 TriFace #4: Points(1 1 2) Neighbors(2 1 3)
6 DualEdges:
 DualEdge 1 : (Δ3)-(1)-------(3)-(Δ2)
 DualEdge 2 : (Δ1)-(1)-------(2)-(Δ2)
 DualEdge 3 : (Δ2)-(1)-------(1)-(Δ4)
 DualEdge 4 : (Δ4)-(2)-------(3)-(Δ1)
 DualEdge 5 : (Δ1)-(2)-------(3)-(Δ3)
 DualEdge 6 : (Δ3)-(2)-------(3)-(Δ4)
```

Check if the 4th edge (DualEdge 4) can be flipped:

```jldoctest DDD
julia> is_flippable(D, 4)
true
```

Flip said edge:

```jldoctest DDD
julia> flip!(D, 4)
DeltaComplex on orientable surface of genus 1 with 2 points
4 TriFaces:
 TriFace #1: Points(1 2 1) Neighbors(3 3 4)
 TriFace #2: Points(1 1 1) Neighbors(4 4 3)
 TriFace #3: Points(1 1 2) Neighbors(2 1 1)
 TriFace #4: Points(1 1 1) Neighbors(2 1 2)
6 DualEdges:
 DualEdge 1 : (Δ3)-(1)-------(3)-(Δ2)
 DualEdge 2 : (Δ4)-(1)-------(2)-(Δ2)
 DualEdge 3 : (Δ2)-(1)-------(3)-(Δ4)
 DualEdge 4 : (Δ4)-(2)-------(3)-(Δ1)
 DualEdge 5 : (Δ1)-(1)-------(3)-(Δ3)
 DualEdge 6 : (Δ3)-(2)-------(2)-(Δ1)
```

Randomly flip edges in `D` until the diameter stabilizes:

```julia-repl
julia> randomize!(D)
10300000
julia> D
DeltaComplex on orientable surface of genus 1 with 2 points
4 TriFaces:
 TriFace #1: Points(1 1 1) Neighbors(3 2 2)
 TriFace #2: Points(1 1 1) Neighbors(1 3 1)
 TriFace #3: Points(1 1 1) Neighbors(1 2 4)
 TriFace #4: Points(2 1 1) Neighbors(4 3 4)
6 DualEdges:
 DualEdge 1 : (Δ3)-(3)-------(2)-(Δ4)
 DualEdge 2 : (Δ3)-(2)-------(2)-(Δ2)
 DualEdge 3 : (Δ4)-(3)-------(1)-(Δ4)
 DualEdge 4 : (Δ1)-(2)-------(3)-(Δ2)
 DualEdge 5 : (Δ1)-(1)-------(1)-(Δ3)
 DualEdge 6 : (Δ2)-(1)-------(3)-(Δ1)
```


## Modular `FlipGraph`

Construct the modular flip graph of a torus with 2 labeled points on it:

```jldoctest
julia> G = flipgraph_modular(1,2)
modular FlipGraph with 9 vertices and 8 edges
```

Construct the modular flip graph of a torus with 2 unlabeled points on it:

```jldoctest
julia> G = flipgraph_modular(1,2;labeled_points=false)
modular FlipGraph with 5 vertices and 4 edges
```