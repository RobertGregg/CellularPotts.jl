####################################################
# Determine cell subgraph
####################################################

using LightGraphs
using CellularPotts
M = ModelParameters()
CPM = CellPotts(M)

g = CPM.graph.network
cellIdx = findall(isequal(1), CPM.graph.σ)

#This function attempts to find a subgraph of the a faster, however performance is about the same 😖
function induced_subgraph_mod(g::T, vlist::AbstractVector{U}) where T <: AbstractGraph where U <: Integer

    gSubRaw = view(g.fadjlist, vlist)
    fadjlist = filter.(x -> x ∈ vlist, gSubRaw)
    #d = Dict(vlist .=> eachindex(vlist))

    for list in fadjlist
        for (j,vertex) in enumerate(list)
            list[j] = searchsortedfirst(vlist,vertex) #assumes vlist is sorted
            #list[j] = d[vertex]
        end
    end

    return SimpleGraph( sum(length,fadjlist) ÷ 2, fadjlist)
end

induced_subgraph_mod(g,cellIdx)


####################################################
# Use a macro to create custom cell types
####################################################

#Cells need rules and properties!

#This is from Agents.jl to create a custom struct with some defaults
macro agent(name, base, fields)
    base_type = Core.eval(@__MODULE__, base)
    base_fieldnames = fieldnames(base_type)
    base_types = [t for t in base_type.types]
    base_fields = [:($f::$T) for (f, T) in zip(base_fieldnames, base_types)]
    res = :(mutable struct $(esc(name)) <: AbstractAgent end)
    push!(res.args[end].args, base_fields...)
    push!(res.args[end].args, map(esc, fields.args)...)
    return res
end