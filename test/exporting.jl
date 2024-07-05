@testset "exporting" begin
    D = deltacomplex(2)
    G = flipgraph_modular(D)
    vats = [Dict("Vertex Attribute" => degree(G, i)) for i in 1:nv(G)]
    export_gml("flipG_test.gml", G, add_diameter=true, vertex_attributes=vats)
    @test isfile("flipG_test.gml")==true
    rm("flipG_test.gml")

    G2 = flipgraph_planar(8)

    vats = [Dict("Vertex Attribute" => degree(G2, i)) for i in 1:nv(G2)]
    eats = [Dict("Edge Attribute" => i) for i in 1:ne(G2)]

    export_gml("flipG2_test", G2; add_diameter=true, vertex_attributes=vats, edge_attributes=eats)
    @test isfile("flipG2_test.gml")==true
    rm("flipG2_test.gml")
end