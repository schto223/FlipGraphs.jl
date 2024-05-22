using Random


@testset "Flip Graph" begin
    @testset "mcKay" begin
        for g in 1:10, p in 1:7
            HD = holey_delta_complex(g,p)
            sps_points = mcKay_points(HD)
            @test length(sps_points) >= 1
            @test all(sp -> length(sp)==p , sps_points)
            @test all(sp -> length(unique(sp))==length(sp) , sps_points)

            sps_trifaces = mcKay_vertices(HD)
            @test length(sps_trifaces) >= 1
            @test all(sp -> length(sp)==nv(HD) , sps_trifaces)
            @test all(sp -> length(unique(sp))==length(sp) , sps_trifaces)

            sps_edges = mcKay_edges(HD)
            @test length(sps_edges) >= 1
            @test all(sp -> length(sp)==ne(HD) , sps_edges)
            @test all(sp -> length(unique(sp))==length(sp) , sps_edges)

            @test is_isomorph(HD, HD, modular=true) == true
            @test is_isomorph(HD, HD, modular=false) == true
        end
    end

    @testset "is_isomorph" begin
        HD = holey_delta_complex(5,7)
        HD2 = holey_delta_complex(5,7)
        flip!(HD2, 5)
        @test is_isomorph(HD, HD2) == false
        flip!(HD2, 5)
        @test is_isomorph(HD, HD2) == true

        HD = holey_delta_complex(8,6)
        HD2 = holey_delta_complex(8,6)
        Random.seed!(71924)
        a = rand(1:ne(HD), 100)
        dir = rand(Bool, 100)
        for i in eachindex(a)
            flip!(HD, a[i], left=dir[i])
        end

        for i in reverse(eachindex(a))
            flip!(HD, a[i], left=!dir[i])
        end
        @test is_isomorph(HD, HD2)
    end

    @testset "flipGraph" begin
        HD = holey_delta_complex(1)
        G = flip_graph(HD, 6, modular=false)
        @test nv(G) >= 11
        @test ne(G) >= 10

        HD = holey_delta_complex(1,2)
        G1 = flip_graph(HD, 2, modular=false)
        G2 = flip_graph(HD, 2, modular=true)
        @test diameter(G1) <= 4
        @test nv(G1) >= nv(G2)
        @test ne(G1) >= ne(G2)
    end

end