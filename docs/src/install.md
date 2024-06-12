# Installation

```@meta
DocTestSetup = quote
    using FlipGraphs
end
```

To install the package, you first need to install [Julia](https://julialang.org). 
After starting Julia, type the following:

```julia-repl
julia> using Pkg 
julia> Pkg.add("FlipGraphs")
```

You can start using the package as follows:

```julia-repl
julia> using FlipGraphs
julia> D = deltacomplex(4,10);
DeltaComplex on orientable surface of genus 4 with 10 points
32 TriFaces:
 TriFace #1: Points(1 1 2) Neighbors(5 15 16)
 TriFace #2: Points(1 1 3) Neighbors(20 17 18)
 ⋮
 TriFace #31: Points(1 1 10) Neighbors(10 32 9)
 TriFace #32: Points(1 1 10) Neighbors(8 9 31)
48 DualEdges:
 DualEdge 1 : (Δ15)-(1)-------(1)-(Δ18)
 DualEdge 2 : (Δ2)-(1)-------(1)-(Δ20)
 ⋮
 DualEdge 47 : (Δ9)-(2)-------(3)-(Δ31)
 DualEdge 48 : (Δ31)-(2)-------(3)-(Δ32)
```

If you need help understanding what a function does or what a structure represents, you can put a question mark in front of it:
```@jldoctest
julia> ?FlipGraph
  struct FlipGraph <: AbstractGraph{Int}

  A Graph representing the flipgraph of a Δ-Complex.

  Vertices are different triangulations of the same surface.
  Two vertices are linked by an edge, if the respective graphs differ only by a single flip.
```
