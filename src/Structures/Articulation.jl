####################################################
# Utility to get articulation points
####################################################

mutable struct Articulation{T<:Integer}
    queue::Vector{T}
    visited::BitVector

    function Articulation(n::T) where T<:Integer

        queue = T[]
        sizehint!(queue, n)

        visited = falses(n)

        return new{T}(queue,visited)
    end
end

function resetArticulation!(F::Articulation)
    
    empty!(F.queue) 
    F.visited .= false

    return nothing
end