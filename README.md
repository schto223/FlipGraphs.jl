# FlipGraphs
*A julia programm for triangulations on closed surfaces and their resulting flip-graphs.*

[pkgeval-img]: https://juliaci.github.io/NanosoldierReports/pkgeval_badges/E/Example.svg
[pkgeval-url]: https://juliaci.github.io/NanosoldierReports/pkgeval_badges/E/Example.html
[docs-img-stable]: https://img.shields.io/badge/docs-stable-green.svg
[docs-img-dev]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-url-stable]: https://schto223.github.io/FlipGraphs.jl/stable
[docs-url-dev]: https://schto223.github.io/FlipGraphs.jl/dev
[codecov-img]:https://codecov.io/gh/schto223/FlipGraphs.jl/graph/badge.svg?token=O2U52538TU
[codecov-url]:https://codecov.io/gh/schto223/FlipGraphs.jl

| **Documentation**                                                         | **Build Status**                                      |
|:-------------------------------------------------------------------------:|:-----------------------------------------------------:|
| [![][docs-img-stable]][docs-url-stable]  [![][docs-img-dev]][docs-url-dev] | [![PkgEval][pkgeval-img]][pkgeval-url] [![codecov][codecov-img]][codecov-url] |


This package is part of my master's thesis on FlipGraphs. It is still unfinished and in active developpment.
If you are not yet familiar with what flip graphs are, I recommend you have a look at [Wikipedia](https://en.wikipedia.org/wiki/Flip_graph), or the documentation for this package.

## Installation

To install the package, you first need to install [Julia](https://julialang.org). 
After starting Julia, type the following:

```julia-repl
julia> using Pkg 
julia> Pkg.add("FlipGraphs")
```

You can start using the package as follows:

```julia-repl
julia> using FlipGraphs
julia> flipgraph_planar(13)
FlipGraphPlanar with 58786 vertices and 293930 edges
```

If you need help with a function, you can put a question mark in front of it,or have a look at the [documentation](https://schto223.github.io/FlipGraphs.jl/dev). 
