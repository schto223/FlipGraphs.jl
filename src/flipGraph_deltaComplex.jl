import Base.LinkedList, Base.LinkedListItem

struct Crossing
    hole_id::Int
    edge_id::Int
    key_holeposition::Int #the key to get the 
    going_in::Bool #true if the edge goes into the hole, false if it comes out
end

struct Hole
    id :: Int
    crossings::LinkedList{Crossing}
    crossingDict::Dict{Int, LinkedListItem}
end

struct holeyDeltaComplex
    D::DeltaComplex
    holes::Vector{Hole}
    edgeCrossings::Vector{Vector{Crossing}}

end


function holeyDeltaComplex(D::DeltaComplex)
    holes = Vector{Hole}(undef, genus(D))
    edgeCrossings = Vector{Vector{Crossing}}([[] for i in 1:ne(D)])
    #TODO this currently assumes that there are at least genus(D) edges where both endpoints are identical
    loopEdges = filter(e->points(e)[1]=points(e)[2] ,edges(D)) 
    for i in 1:length(holes) 
        loop = loopEdges[i]
        T_start = get_vertex(D,vertices(e)[1])
        while true
            e = left_edge(T_start, loop)
            #c = Crossing(i, e)
        end
    end

end