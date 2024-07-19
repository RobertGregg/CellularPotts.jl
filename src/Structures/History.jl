####################################################
# Logging Function
####################################################


struct History{T<:Integer, C, N}
    space::CellSpace{T,C,N}
    counter::Vector{T}      #Time step where space was changed
    node::Vector{T}         #location in space that was changed
    id::Vector{T}           #new id for that location
    type::Vector{T}         #new type for that location
end

History(space) = History(space,Int[],Int[],Int[],Int[])