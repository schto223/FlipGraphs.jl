# Installation

To install the package, you first need to install [Julia](https://julialang.org). 
After starting Julia, type the following:

```julia-repl
julia> using Pkg 
julia> Pkg.add("FlipGraphs")
```

You can start using the package as follows:

```julia-repl
julia> using FlipGraphs
julia> createDeltaComplex(4,10)
```

You can get help for a function or structure by putting a question mark in front.
```julia-repl
help?> FlipGraph
  struct FlipGraph <: AbstractGraph{Int}

  A Graph representing the flipgraph of a Δ-Complex.

  Vertices are different triangulations of the same surface.
  Two vertices are linked by an edge, if the respective graphs differ only by a single flip.
```
