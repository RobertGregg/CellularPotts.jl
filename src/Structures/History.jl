####################################################
# Logging Function
####################################################

#TODO Why does History need penalty?
struct History{T<:Integer, C, N}
    space::CellSpace{T,C,N}
    counter::Vector{T}
    position::Vector{T}
    nodeID::Vector{T}
    nodeType::Vector{T}
end

History(space) = History(space,Int[],Int[],Int[],Int[])