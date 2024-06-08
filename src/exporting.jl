
export export_gml

"""
    export_gml(fpn::String, G::AbstractGraph{::Integer})

Save the graph `G` as a .gml file.

`G` can be either a `FlipGraph` or `FlipGraphPlanar`.

# Example
```julia-repl
julia> G = flipgraph_planar(10);
julia> export_gml("C:/Users/USERNAME/Desktop/filename.gml", G);
```
"""
function export_gml(fpn::String, G::AbstractGraph{<:Integer})
    open(fpn,"w") do file
        write(file,"graph [")
        write(file,"\n\tdirected 0")
        for i in eachindex(G.V)
            write(file,"\n\tnode [")
            write(file,"\n\t\tid $i")
            write(file,"\n\t\tlabel \"$i\"")
            write(file,"\n\t]")
        end
        for e in edges(G)
            write(file,"\n\tedge [")
            write(file,"\n\t\tsource ", string(e.src))
            write(file,"\n\t\ttarget ", string(e.dst))
            write(file,"\n\t]")
        end
        write(file,"\n]")
    end
end