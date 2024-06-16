@testset "show" begin
    io = IOBuffer()

    @testset "convex polygon" begin
        p = triangulated_polygon(10)
        show(io,"text/plain",p)
        p = triangulated_polygon(30)
        show(io,"text/plain",p)

        G = flipgraph_planar(10)
        show(io,"text/plain",G)
    end

    @testset "closed surfaces" begin
        D = deltacomplex(3,4)
        show(io,"text/plain",D)
        show(io,"text/plain",D.V[1])
        show(io,"text/plain",D.E[2])

        D = deltacomplex_non_orientable(5,2)
        show(io,"text/plain",D)
        show(io,"text/plain",D.V[1])
        show(io,"text/plain",D.E[2])

        D = deltacomplex(1)
        G = flipgraph_modular(D)
        show(io,"text/plain",G)
    end

    
end