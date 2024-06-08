@testset "exporting" begin
    D = deltacomplex(1)
    G = flipgraph_modular(D, 10)

    export_gml("flipG_test", G)
    @test isfile("flipG_test")==true
    rm("flipG_test")
end