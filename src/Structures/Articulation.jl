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

        return new(queue,visited)
    end
end

function resetArticulation!(F::Articulation)
    
    empty!(F.queue) 
    F.visited .= false

    return nothing
end

#TODO Remove out of this file
#Returns true is deleting the target node would disconnect the graph 
function isfragmented(cpm)

    #Unpack cpm and reset state (no allocations right?)
    g = cpm.space
    step = cpm.step
    F = cpm.fragment

    #TODO Not checking Medium mean cells can enclose bits of medium
    #No need to check Medium for articulation points  
    if iszero(step.targetCellID)
        return false
    end

    if g.neighborCount == 4
        return isfragmented_quick(g, step.targetNode, step.targetCellID)
    end
    
    resetArticulation!(F)
    
    #find a neighbor to the target node that has the same ID as the target
    neighborTargetCellID = findfirst(isequal(step.targetCellID), g.nodeIDs[i] for i in step.targetNeighborNodes)
    neighborTargetNode = step.targetNeighborNodes[neighborTargetCellID]

    #Mark the target neighbor as visited
    F.visited[neighborTargetNode] = true

    #Add the target neighbor to the queue
    push!(F.queue, neighborTargetNode)

    #Breath first search through the cell
    while !isempty(F.queue)
        currentNodeID = popfirst!(F.queue)

        #Loop through neighbors in target cell and not equal targetCellID
        for neighbor in Iterators.filter(n -> n ≠ step.targetNode && g.nodeIDs[n] == step.targetCellID, neighbors(g,currentNodeID))
            if !F.visited[neighbor]
                F.visited[neighbor] = true
                push!(F.queue, neighbor)
            end

            # #TODO Need to check this, but if we reach all the neighbors of targetNode then I think we can stop
            if all(F.visited[n] for n in step.targetNeighborNodes)
                return false
            end
        end
    end

    #TODO Use state volume? This won't work for Medium.
    #If we visited all the nodes (except targetNode) then cell is still connected
    if count(isequal(step.targetCellID), g.nodeIDs) == count(F.visited) + 1
        return false
    end

    return true
end

#Special case for 2D, planar spaces
#Count the number of vertices, edges, and faces changed when removing "node" from the graph
#V-E+F == 2 for all planar connected graphs ∴ we should see no change when removing "node"
function isfragmented_quick(g, node, cellID)
    
    V = 1
    
    E = count(g.nodeIDs[n] == cellID for n in neighbors(g,node))
    
    F = 0
    for (x,y) in allpairs(neighbors(g,node))
        if count(n -> insorted(n, neighbors(g,x)), neighbors(g,y)) == 2
            F += 1
        end
    end

    if V-E+F == 0
        return false
    else
        return true
    end
end
    
allpairs(v) = Iterators.filter(i -> isless(i...), Iterators.product(v,v))



# function testfragmented(g, node)
    
#     V = 1
    
#     E = length(neighbors(g,node))
    
#     F = 0
#     for (x,y) in allpairs(neighbors(g,node))

#         # O(n²) but still faster for small n
#         if count(n -> any(isequal(n), neighbors(g,x)), neighbors(g,y)) == 2
#             F += 1
#         end
#     end

#     if V-E+F == 0
#         return "not articulation"
#     else
#         return "articulation"
#     end
# end

# allpairs(v) = Iterators.filter(i -> isless(i...), Iterators.product(v,v))

# using Graphs
# # 4 and 7 are articulation points
# g = SimpleGraph(8)
# add_edge!(g,1,2)
# add_edge!(g,1,3)
# add_edge!(g,2,4)
# add_edge!(g,3,4)
# add_edge!(g,4,5)
# add_edge!(g,4,6)
# add_edge!(g,5,7)
# add_edge!(g,6,7)
# add_edge!(g,7,8)


# for v in vertices(g)

#     out = testfragmented(g, v)

#     println("$(v) → $(out) point")
# end