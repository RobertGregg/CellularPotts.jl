####################################################
# Determine cell subgraph
####################################################

using Graphs
using CellularPotts

g = CellSpace(50,50)
cellIdx = rand(1:50^2,100) |> sort

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

g1 = induced_subgraph_mod(g,cellIdx)



function DFS(g::SimpleGraph)
    src = first(vertices(g))
    n = nv(g)
    stack = Int[]
    visited = falses(n)
    depth = zeros(Int,n)
    low = zeros(Int,n)

    push!(stack, src)

    while !isempty(stack)
        currentNode = pop!(stack)
        if !visited[currentNode]
            visited[currentNode] = true
            for neighbor in neighbors(g, currentNode)
                push!(stack, neighbor)
            end
        end
    end
end
