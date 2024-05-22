using Documenter, FlipGraphs#, DocumenterInterLinks
import Graphs: Edge

#links = InterLinks(
#    "Graphs" => "https://juliagraphs.org/Graphs.jl/stable/",
#    "Julia" => (
#        "https://docs.julialang.org/en/v1/",
#        joinpath(@__DIR__, "inventories", "Julia.toml")
#    ),
#    "Documenter" => (
#        "https://documenter.juliadocs.org/stable/",
#        "https://documenter.juliadocs.org/stable/objects.inv",
#        joinpath(@__DIR__, "inventories", "Documenter.toml")
#    )
#);

makedocs(
    sitename="FlipGraphs Documentation", 
    modules= [FlipGraphs],
    doctest = false, #remove later
    pages = [
        "About" => "index.md",
        "Installation" => "install.md",
        "Quick Start" => "quickStart.md",
        "Convex Polygons" => ["polygonTriangulation.md", "flipGraph_planar.md"],
        "Closed Surfaces" => ["deltaComplex.md", "holeyDeltaComplex.md", "flipGraph.md"],
        "General Utilities" => "generalUtilities.md"
        #"Plotting" => ["plotting.md"]
    ],
    checkdocs = :none #:exports 
)

deploydocs(
    repo = "github.com/schto223/FlipGraphs.jl.git",
)
