####################################################
# Logging Function
####################################################

#TODO Why does History need space?
struct History{T<:Integer, C, N}
    space::CellSpace{T,C,N}
    step::Vector{T}
    idx::Vector{T}
    nodeID::Vector{T}
    nodeType::Vector{T}
    penalty::Vector{Vector{Float64}}
end

History(space) = History(space,Int[],Int[],Int[],Int[], Vector{Float64}[])