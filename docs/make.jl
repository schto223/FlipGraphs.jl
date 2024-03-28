using Documenter, FlipGraphs

#makedocs(sitename="FlipGraphs Documentation", modules= [FlipGraphs])

makedocs(
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