@testset "exporting" begin
    D = deltacomplex(2)
    G = flipgraph_modular(D)
    export_gml("flipG_test.gml", G, :diameter)
    @test isfile("flipG_test.gml")==true
    rm("flipG_test.gml")

    G2 = flipgraph_planar(8)
    export_gml("flipG2_test", G2, :diameter)
    @test isfile("flipG2_test.gml")==true
    rm("flipG2_test.gml")
end