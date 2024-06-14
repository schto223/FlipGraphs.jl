using Random
import Graphs: SimpleEdge, Edge

function is_isomorphic_safe(D1::DeltaComplex, D2::DeltaComplex; fix_points::Bool = true)
    D=deepcopy(D1)
    if !fix_points
        rename_points!(D,mcKay_points(D)[1])
        point_perms = mcKay_points(D2)
    else
        point_perms = [collect(1:np(D2))]
    end

    rename_vertices!(D,mcKay_vertices(D)[1])
    rename_edges!(D,mcKay_edges(D)[1])

    p = 1
    for pp in point_perms
        Dp = rename_points!(deepcopy(D2), pp)
        v=1
        for pv in mcKay_vertices(Dp)
            Dv = rename_vertices!(deepcopy(Dp), pv)
            e=1
            for pe in mcKay_edges(Dv)
                De = rename_edges!(deepcopy(Dv), pe)
                bo = true
                for i in 1:nv(D)
                    if !is_equivalent(D.V[i], De.V[i])
                        bo=false
                        break
                    end
                end
                if bo
                    return true
                end
                e+=1
            end
            v+=1
        end
        p+=1
    end
    return false
end


@testset "Flip Graph" begin
    @testset "mcKay" begin
        for g in 1:9, p in 1:7
            D = deltacomplex(g,p)
            ps_points = mcKay_points(D)
            @test length(mcKay_points(D, only_one=true)) == 1
            @test mcKay_points(D, only_one=true)[1] in ps_points
            @test length(ps_points) >= 1
            @test all(sp -> length(sp)==p , ps_points)
            @test all(sp -> length(unique(sp))==length(sp) , ps_points)

            ps_vertices = mcKay_vertices(D)
            @test length(mcKay_vertices(D, only_one=true)) == 1
            @test mcKay_vertices(D, only_one=true)[1] in ps_vertices
            @test length(ps_vertices) >= 1
            @test all(sp -> length(sp)==nv(D) , ps_vertices)
            @test all(sp -> length(unique(sp))==length(sp) , ps_vertices)

            ps_edges = mcKay_edges(D)
            @test length(mcKay_edges(D, only_one=true)) == 1
            @test mcKay_edges(D, only_one=true)[1] in ps_edges
            @test length(ps_edges) >= 1
            @test all(sp -> length(sp)==ne(D) , ps_edges)
            @test all(sp -> length(unique(sp))==length(sp) , ps_edges)

            @test is_isomorphic(D, D, labeled_points=false) == true
            @test is_isomorphic(D, D, labeled_points=true) == true
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

    @testset "FlipGraph" begin
        D = deltacomplex(1)
        G = flipgraph_modular(D)
        @test nv(G) == 1
        @test ne(G) == 0

        D = deltacomplex(1,2)
        G1 = flipgraph_modular(D, labeled_points=true)
        G2 = flipgraph_modular(D, labeled_points=false)
        @test diameter(G1) == 8
        @test nv(G1) >= nv(G2)
        @test ne(G1) >= ne(G2)

        G = flipgraph_modular(1,2)
        @test edgetype(G) == SimpleEdge{Int32}
        @test is_directed(G) == false
        @test is_directed(FlipGraph) == false
        @test has_edge(G, Edge(1,2)) == true
        @test has_edge(G, get_vertex(G,1), get_vertex(G,2)) == true
        @test 1 in outneighbors(G, 2)
        @test 2 in inneighbors(G, 1)
        V = vertices(G)
        Vd = vertices_deltacomplex(G)
        for i in eachindex(V)
            @test V[i].D == Vd[i]
            for j in eachindex(V)
                if i!=j
                    @test is_isomorphic(Vd[i], Vd[j]) == false
                end
            end
        end

    end

    @testset "flipgraph_modular" begin
        nvs_labeled = [[1,9,236],[9,713]]
        nes_labeled = [[0,8,591],[13,2938]]
        nvs = [[1,5,46,669],[9,368]]
        nes = [[0,4,98,2684],[13,1478]]
        #labeled
        for g in eachindex(nvs_labeled)
            for p in eachindex(nvs_labeled[g])
                G = flipgraph_modular(g, p, labeled_points=true)
                @test nv(G) == nvs_labeled[g][p]
                @test ne(G) == nes_labeled[g][p]
            end
        end

        #unlabeled
        for g in eachindex(nvs)
            for p in eachindex(nvs[g])
                G = flipgraph_modular(g, p, labeled_points=false)
                @test nv(G) == nvs[g][p]
                @test ne(G) == nes[g][p]
            end
        end

        #sphere
        for p in 3:5
            G = flipgraph_modular(0, p, labeled_points=true)
            @test nv(G) == [4, 64, 2240][p-2]
            @test ne(G) == [3, 120, 7200][p-2]
        end
        for p in 3:7
            G = flipgraph_modular(0, p, labeled_points=false)
            @test nv(G) == [2, 6, 26, 191, 1904][p-2]
            @test ne(G) == [1, 5, 48, 658, 9627][p-2]
        end
    end
end


function isGood(fgvc::FGVertexCandidate)
    D = fgvc.D
    p = 1
    for pp in fgvc.point_perms
        Dp = rename_points!(deepcopy(D), pp)
        if !sameDcomplex(Dp, get_D(fgvc, p))
            println("badP")
        end
        v=1
        for pv in get_vertex_perms(fgvc, p)
            Dv = rename_vertices!(deepcopy(Dp), pv)
            if !sameDcomplex(Dv, get_D(fgvc, p, v))
                println("badV")
            end
            e=1
            for pe in get_edge_perms(fgvc, p, v)
                De = rename_edges!(deepcopy(Dv), pe)
                if !sameDcomplex(De, get_D(fgvc, p, v, e))
                    println("badE")
                end
                e+=1
            end
            v+=1
        end
        p+=1
    end
end

