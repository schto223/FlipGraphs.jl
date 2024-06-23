"""
    export_gml(fpn::String, G::AbstractGraph{::Integer}, kwargs...)

Save the graph `G` as a .gml file.

`G` can be either a `FlipGraph` or `FlipGraphPlanar`.

# Examples
```julia-repl
julia> G = flipgraph_planar(10);
julia> export_gml("C:/Users/USERNAME/Desktop/filename.gml", G);
```

By adding the additional symbol `:diameter`, the nodes get a value *diameter* which corresponds to the diameter of the `DeltaComplex` or `TriangulatedPolygon` it models.
Be aware however, that this diameter is computed on the run and will therefore significantly slow donw this export method.

```julia-repl
julia> G = flipgraph_modular(1,3,labeled_points = true);
julia> export_gml("C:/Users/USERNAME/Desktop/filename.gml", G, :diameter);
```

"""
function export_gml(fpn::String, G::AbstractGraph, kwargs...)
    if '.' in fpn
        fpn[end-2:end] == "gml" || throw(ArgumentError("The filename should end in \".gml\""))
    else
        fpn *= ".gml"
    end
    open(fpn,"w") do file
        write(file,"graph [")
        write(file,"\n\tdirected 0")
        if typeof(G) == FlipGraph
            write(file,"\n\tlabeled_points $(G.labeled_points)")
        elseif typeof(G) == FlipGraphPlanar
            write(file,"\n\tmodular $(G.modular)")
        end
        for i in eachindex(vertices(G))
            write(file,"\n\tnode [")
            write(file,"\n\t\tid $i")
            write(file,"\n\t\tlabel \"$i\"")
            if :diameter in kwargs
                write(file,"\n\t\tdiameter $(diameter(G.V[i]))")
            end
            write(file,"\n\t]")
        end
        for e in edges(G)
            write(file,"\n\tedge [")
            write(file,"\n\t\tsource ", string(e.src))
            write(file,"\n\t\ttarget ", string(e.dst))
            write(file,"\n\t]")
        end
        write(file,"\n]");
    end
    open(io->read(io,String), fpn);
end