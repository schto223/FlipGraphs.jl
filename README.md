# FlipGraphs
*A julia programm for triangulations on closed surfaces and their resulting flip-graphs.*
<!--
| **Documentation**                                                         | **Build Status**                                      |
|:-------------------------------------------------------------------------:|:-----------------------------------------------------:|
| [![][docs-stable-img]][docs-stable-url] [![][docs-dev-img]][docs-dev-url] | [![][ga-img]][ga-url] [![][codecov-img]][codecov-url] |
-->
<!--[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://schto223.github.io/FlipGraphs.jl/stable)-->
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://schto223.github.io/FlipGraphs.jl/dev)
[![codecov](https://codecov.io/gh/schto223/FlipGraphs.jl/graph/badge.svg?token=O2U52538TU)](https://codecov.io/gh/schto223/FlipGraphs.jl)

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
julia> createDeltaComplex(4,10)
```

You can get help for a function by putting a question mark in front.
