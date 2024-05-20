using Documenter, FlipGraphs, Graphs

makedocs(
    sitename="FlipGraphs Documentation", 
    modules= [FlipGraphs],
    doctest = false, #remove later
    pages = [
        "About" => "index.md",
        "Convex Polygons" => ["polygonTriangulation.md", "flipGraph_planar.md"],
        "Closed Surfaces" => ["deltaComplex.md", "holeyDeltaComplex.md", "flipGraph.md"]#,
        #"Plotting" => ["plotting.md"]
    ],
    checkdocs = :exports # remove later or replace by :exports
)

"""makedocs(
    modules = [FlipGraphs],
    format = Documenter.HTML(
        # Use clean URLs, unless built as a "local" build
        prettyurls = !("local" in ARGS),
        #canonical = "https://juliadocs.github.io/Documenter.jl/stable/",
    ),
    sitename = "FlipGraphs Documentation",
#    authors = "Tom Schmit",
    pages = [
        "About" => "index.md"
    ]
)

deploydocs(
    repo = "github.com/schto223/FlipGraphs.jl.git",
)
"""