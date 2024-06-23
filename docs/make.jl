using Documenter, FlipGraphs#, DocumenterInterLinks
import Graphs: Edge

#links = InterLinks(
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

#DocMeta.setdocmeta!(FlipGraphs, :DocTestSetup, :(using FlipGraphs); recursive=true)

makedocs(
    sitename="FlipGraphs.jl",     
    modules= [FlipGraphs],
    #format = Documenter.LaTeX(platform = "none"), #format = Documenter.LaTeX(),
    #format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    doctest = true,
    pages = [
        "Home" => ["index.md", "install.md", "quickStart.md"],
        "Convex Polygons" => ["polygonTriangulation.md", "flipGraphPlanar.md"],
        "Closed Surfaces" => ["deltaComplex.md", "flipGraph.md"],
        "General Utilities" => ["generalUtilities.md", "exporting.md"]
        #"Plotting" => ["plotting.md"]
    ],
    #checkdocs = :exports 
    checkdocs = :none
)

deploydocs(;
    repo = "github.com/schto223/FlipGraphs.jl.git"#,
    #devbranch = "dev",
    #versions = "v^"
)
