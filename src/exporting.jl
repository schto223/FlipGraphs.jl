"""
    export_gml(fpn::String, G::AbstractGraph{::Integer}, kwargs...)

Save the graph `G` as a .gml file.

`G` can be either a `FlipGraph` or `FlipGraphPlanar`.

# Arguments
- vertex_attributes::Vector{Dict} : If provided, the `i`-th vertex gets the attributes from the `i`-th `Dict`. 
                                    Each `key`, `value` pair defines an attribute, where the `key` is the name of the attribute and `value` is its value.
- edge_attributes::Vector{Dict}   : If provided, the `i`-th edge gets the attributes from the `i`-th `Dict`. 
                                    Each `key`, `value` pair defines an attribute, where the `key` is the name of the attribute and `value` is its value.
- add_diameter = false :            If set to true, every vertex in the exported tree gets an attribute called diameter with the respective diameter of the `DeltaComplex` of that vertex.

# Examples
```julia-repl
julia> G = flipgraph_planar(10);
julia> export_gml("C:/Users/USERNAME/Desktop/filename.gml", G);
```

By adding the additional symbol `:diameter`, the nodes get a value *diameter* which corresponds to the diameter of the `DeltaComplex` or `TriangulatedPolygon` it models.
Be aware however, that this diameter is computed on the run and will therefore significantly slow down this export method.

```julia-repl
julia> G = flipgraph_modular(1,3,labeled_points = true);
julia> export_gml("C:/Users/USERNAME/Desktop/filename.gml", G, add_diameter=true);
```

"""
function export_gml(fpn::String, G::AbstractGraph; vertex_attributes::Vector=[], edge_attributes::Vector=[], add_diameter = false)
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
            if add_diameter
                write(file,"\n\t\tdiameter $(diameter(G.V[i]))")
            end
            if !isempty(vertex_attributes)
                for (key,value) in vertex_attributes[i]
                    write(file,"\n\t\t$(key) $(value)")
                end
            end
            write(file,"\n\t]")
        end
        i = 0
        for e in edges(G)
            i += 1
            write(file,"\n\tedge [")
            write(file,"\n\t\tsource ", string(e.src))
            write(file,"\n\t\ttarget ", string(e.dst))
            if !isempty(edge_attributes)
                for (key,value) in edge_attributes[i]
                    write(file,"\n\t\t$(key) $(value)")
                end
            end
            write(file,"\n\t]")
        end
        write(file,"\n]");
    end
    open(io->read(io,String), fpn);
end