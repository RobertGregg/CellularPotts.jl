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
# Determine if a node is an Articulation point
####################################################

#Currently, articulation points are recalculated each time MHStep is called which is very slow/wasteful
#Need a better way to calculate to determine how articulation points change after a node is added removed