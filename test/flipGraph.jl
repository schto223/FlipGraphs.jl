@testset "Flip Graph" begin
    @testset "mcKay" begin
        for g in 1:10, p in 1:10
            HD = holeyDeltaComplex(g,p)
            sps_points = mcKay_points(HD)
            @test length(sps_points) >= 1
            @test all(sp -> length(sp)==p , sps_points)
            @test all(sp -> length(unique(sp))==length(sp) , sps_points)

            sps_trifaces = mcKay_triFaces(HD)
            @test length(sps_trifaces) >= 1
            @test all(sp -> length(sp)==nv(HD) , sps_trifaces)
            @test all(sp -> length(unique(sp))==length(sp) , sps_trifaces)

            sps_edges = mcKay_edges(HD)
            @test length(sps_edges) >= 1
            @test all(sp -> length(sp)==ne(HD) , sps_edges)
            @test all(sp -> length(unique(sp))==length(sp) , sps_edges)

            @test is_isomorph(HD, HD, false) == true
            @test is_isomorph(HD, HD, true) == true
        end
    end

    @testset "is_isomorph" begin
        HD = holeyDeltaComplex(5,7)
        HD2 = holeyDeltaComplex(5,7)
        flip!(HD2, 5, true)
        @test is_isomorph(HD, HD2) == false
        flip!(HD2, 5, true)
        @test is_isomorph(HD, HD2) == true


        HD = holeyDeltaComplex(9,12)
        HD2 = holeyDeltaComplex(9,12)
        Random.seed!(71924)
        a = rand(1:ne(HD), 1000)
        dir = rand(Bool, 1000)
        for i in eachindex(a)
            flip!(HD, a[i], dir[i])
        end

        for i in reverse(eachindex(a))
            flip!(HD, a[i], !dir[i])
        end
        @test is_isomorph(HD, HD2)
    end

    @testset "flipGraph" begin
        HD = holeyDeltaComplex(1)
        G = construct_FlipGraph(HD, 10 ,true)
        @test nv(G) == 11
        @test ne(G) == 10

        HD = holeyDeltaComplex(5,9)
        G1 = construct_FlipGraph(HD, 10 ,true)
        G2 = construct_FlipGraph(HD, 10 ,false)
        @test diameter(G1) <= 20
        @test nv(G1) >= nv(G2)
        @test ne(G1) >= ne(G2)
    end

end