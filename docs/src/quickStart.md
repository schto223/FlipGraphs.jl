# Quick Start

If you're already familiar with the concept of **flipgraphs**, **triangulations** on **closed surfaces** and **Δ-complexes**, and don't to read the whole documentation of what what is, then here are some quick examples of what you can do with this package.

In any other case, please be sure to have a look at the rest of the documentation first.

## Triangulated Convex Polygon

Create a triangulated convex 10-gon:

```julia-repl
julia> g = triangulated_polygon(10)
TriangulatedPolygon with 10 vertices, and adjacency list:
[[2, 10], [1, 3, 10], [2, 4, 10, 9], [3, 5, 9, 8], [4, 6, 8, 7], [5, 7], [6, 8, 5], [7, 9, 4, 5], [8, 10, 3, 4], [9, 1, 2, 3]]
```

Check if the edge going from vertex 2 to vertex 10 can be flipped:

```julia-repl
julia> is_flippable(g, 2, 10)
true
```

Flip said edge:

```julia-repl
julia> flip!(g, 2, 10)
TriangulatedPolygon with 10 vertices, and adjacency list:
[[2, 10, 3], [1, 3], [2, 4, 10, 9, 1], [3, 5, 9, 8], [4, 6, 8, 7], [5, 7], [6, 8, 5], [7, 9, 4, 5], [8, 10, 3, 4], [9, 1, 3]]
```

Construct the flipgraph of an Octagon:

```julia-repl
julia> G = flipgraph_planar(8)
FlipGraphPlanar with 132 vertices and 330 edges
```


## Δ-Complex / Triangulation of closed surface

### `DeltaComplex`

A `DeltaComplex` is a simplified triangulation on a closed surface.
(Can be used to compute things like the diameter etc. But does not offer a unique modelisation of a triangulation on a closed surface)

Create a `DeltaComplex` of a surface of genus 1 with 2 points 

```julia-repl
julia> D = delta_complex(1, 2)
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

```julia-repl
julia> is_flippable(D, 4)
true
```

Flip said edge:

```julia-repl
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

Randomly flip edges in `D` until it the diameter stabilizes:

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


### `HoleyDeltaComplex`

A `DeltaComplex` is a more elaborate modelisation of a triangulation on a closed surface. This allows to determine if two triangulisations are homeomorph to each other or not. This comes however at the cost of efficiency.

```julia-repl
julia> HD = holey_delta_complex(1, 2)
HoleyDeltaComplex on orientable surface of genus 1 with: 2 Points; 4 TriFaces; 
6 DualEdges; 1 Hole:
Hole 1 : --<--3⤈-1⤈-6⤉--<--
```

Construct a local image of the flipgraph containing `HD` up to a depth of 5.

```julia-repl
julia> flip_graph(HD, 5)
FlipGraph with 98 vertices and 154 edges
```