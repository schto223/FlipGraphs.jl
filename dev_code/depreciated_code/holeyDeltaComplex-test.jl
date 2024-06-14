using Random

function holes_intact(HD::HoleyDeltaComplex)
    for H in HD.holes
        c_first = get_crossing(H)
        c = c_first
        skipfirst = true
        while c != c_first || skipfirst
            skipfirst = false
            if c.previous.next != c || c.next.previous != c
                return false
            end
            d_prev = get_edge(HD, c.previous.edge_id)
            d = get_edge(HD, c.edge_id)
            if !(d.triangles[1] in vertices_id(d_prev) || d.triangles[2] in vertices_id(d_prev))
                return false
            end
            c = c.next
        end
    end
    return true
end  

@testset "HoleyDeltaComplex" begin

    @testset "Sphere" begin
        for p = 3:10
            HD = holey_delta_complex(0, p)
            @test nv(HD) == 2*(p-2)  
            @test ne(HD) == 3*(p-2) 
            @test np(HD) == p
            @test euler_characteristic(HD) == 2
            @test genus(HD) == 0
            @test isempty(HD.holes)
        end
    end

    @testset "Orientable surfaces" begin
        for g = 1:8
            for p = 1:8
                HD = holey_delta_complex(g, p)
                @test genus(HD) == g
                @test np(HD) == p
                @test euler_characteristic(HD) == 2-2*g    
                @test sum(H -> num_crossings(H), HD.holes) == length(reduce(vcat, HD.edge_crossings)) 
                @test holes_intact(HD)   
            end
        end

        #test flip
        HD = holey_delta_complex(3,7)
        @test sum(point_degrees(HD)) == 2*ne(HD)
        d = get_edge(HD,1)
        e_copy = deepcopy(d)
        t1,t2 = vertices_id(d)
        T1 = deepcopy(get_vertex(HD, t1))
        T2 = deepcopy(get_vertex(HD, t2))
        flip!(HD,d)
        @test genus(HD) == 3
        @test np(HD) == 7
        @test euler_characteristic(HD) == -4  
        @test sum(H -> num_crossings(H), HD.holes) == length(reduce(vcat, HD.edge_crossings))  
        @test holes_intact(HD)   
        flip!(HD, d; left=false)        
        @test all(e_copy.triangles.==d.triangles) && all(e_copy.sides==d.sides) && d.is_twisted==e_copy.is_twisted
        @test sum(H -> num_crossings(H), HD.holes) == length(reduce(vcat, HD.edge_crossings))  
        @test holes_intact(HD)   
    end 

   @testset "renaming & is_isomorph" begin
        Random.seed!(123)
        HD = holey_delta_complex(5,5)
        HD2 = deepcopy(HD)
        rename_points!(HD2, shuffle(1:np(HD)))
        @test is_isomorph(HD, HD2; modular=true) == true

        HD2 = deepcopy(HD)
        rename_vertices!(HD2, shuffle(1:nv(HD)))
        @test is_isomorph(HD, HD2) == true

        HD2 = deepcopy(HD)
        rename_edges!(HD2, shuffle(1:ne(HD)))
        @test is_isomorph(HD, HD2) == true
   end

    @testset "diameter & reversability" begin 
        HD = holey_delta_complex(10,20)
        HD2 = holey_delta_complex(10,20)
        @test 1 <= diameter_triangulation(HD) <= np(HD)-1
        @test 1 <= diameter_deltaComplex(HD) <= nv(HD)-1
        Random.seed!(1234)
        a = rand(1:ne(HD), 10)
        dir = rand(Bool, 10)
        for i in eachindex(a)
            flip!(HD, a[i]; left=dir[i])
        end 
        @test 1 <= diameter_triangulation(HD) <= np(HD)-1
        @test 1 <= diameter_deltaComplex(HD) <= nv(HD)-1
        @test sum(H -> num_crossings(H), HD.holes) == length(reduce(vcat, HD.edge_crossings))  
        @test holes_intact(HD)   

        for i in reverse(eachindex(a))
            flip!(HD, a[i]; left=!dir[i])
        end

        @test is_isomorph(HD, HD2) == true
    end
end


  
