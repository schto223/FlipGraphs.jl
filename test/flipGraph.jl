using Random


@testset "Flip Graph" begin
    @testset "mcKay" begin
        for g in 1:10, p in 1:7
            D = deltacomplex(g,p)
            sps_points = mcKay_points(D)
            @test length(sps_points) >= 1
            @test all(sp -> length(sp)==p , sps_points)
            @test all(sp -> length(unique(sp))==length(sp) , sps_points)

            sps_trifaces = mcKay_vertices(D)
            @test length(sps_trifaces) >= 1
            @test all(sp -> length(sp)==nv(D) , sps_trifaces)
            @test all(sp -> length(unique(sp))==length(sp) , sps_trifaces)

            sps_edges = mcKay_edges(D)
            @test length(sps_edges) >= 1
            @test all(sp -> length(sp)==ne(D) , sps_edges)
            @test all(sp -> length(unique(sp))==length(sp) , sps_edges)

            @test is_isomorphic(D, D, fix_points=false) == true
            @test is_isomorphic(D, D, fix_points=true) == true
        end
    end

    @testset "is_isomorphic" begin
        D = deltacomplex(5,7)
        D2 = deltacomplex(5,7)
        flip!(D2, 5)
        @test is_isomorphic(D, D2) == false
        flip!(D2, 5)
        @test is_isomorphic(D, D2) == true

        D = deltacomplex(8,6)
        D2 = deltacomplex(8,6)
        Random.seed!(71924)
        a = rand(1:ne(D), 100)
        dir = rand(Bool, 100)
        for i in eachindex(a)
            flip!(D, a[i], left=dir[i])
        end

        for i in reverse(eachindex(a))
            flip!(D, a[i], left=!dir[i])
        end
        @test is_isomorphic(D, D2)
    end

    @testset "flipGraph" begin
        D = deltacomplex(1)
        G = flipgraph_modular(D, 6)
        @test nv(G) == 2
        @test ne(G) == 1

        D = deltacomplex(1,2)
        G1 = flipgraph_modular(D, 100, fix_points=true)
        G2 = flipgraph_modular(D, 100, fix_points=false)
        @test diameter(G1) == 8
        @test nv(G1) >= nv(G2)
        @test ne(G1) >= ne(G2)
    end

end