@testset "exporting" begin
    D = deltacomplex(2)
    G = flipgraph_modular(D)

    export_gml("flipG_test.gml", G)
    @test isfile("flipG_test.gml")==true
    rm("flipG_test.gml")
end