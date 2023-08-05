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

#TODO Remove out of this file
#Returns true is deleting the target node would disconnect the graph 
function isfragmented(cpm)

    #Unpack cpm and reset state (no allocations right?)
    g = cpm.space
    source = cpm.step.source
    target = cpm.step.target
    F = cpm.fragment

    #TODO Not checking Medium mean cells can enclose bits of medium
    #No need to check Medium for articulation points  
    if iszero(target.id)
        return false
    end

    if maxNeighborCount(cpm.space) == 4
        return isfragmented_quick(g, target.node, target.id)
    end
    
    resetArticulation!(F)
    
    #find a neighbor to the target node that has the same ID as the target
    startID = findfirst(isequal(target.id), g.nodeIDs[i] for i in target.neighbors)
    startNode = target.neighbors[startID]

    #Mark the target neighbor as visited
    F.visited[startNode] = true

    #Add the target neighbor to the queue
    push!(F.queue, startNode)

    #Breath first search through the cell
    while !isempty(F.queue)
        currentNode = popfirst!(F.queue)

        #Loop through neighbors in target cell and not equal target ID
        for neighbor in Iterators.filter(n -> n ≠ target.node && g.nodeIDs[n] == target.id, neighbors(g,currentNode))
            if !F.visited[neighbor]
                F.visited[neighbor] = true
                push!(F.queue, neighbor)
            end

            # #TODO Need to check this, but if we reach all the neighbors of target node then I think we can stop
            if all(F.visited[n] for n in target.neighbors)
                return false
            end
        end
    end

    #TODO Use state volume? This won't work for Medium.
    #If we visited all the nodes (except target node) then cell is still connected
    if count(isequal(target.id), g.nodeIDs) == count(F.visited) + 1
        return false
    end

    return true
end

#Special case for 2D, planar spaces
#Count the number of vertices, edges, and faces changed when removing "node" from the graph
#V-E+F == 2 for all planar connected graphs ∴ we should see no change when removing "node"
function isfragmented_quick(g, node, id)
    
    V = 1
    
    E = count(g.nodeIDs[n] == id for n in neighbors(g,node))
    
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
    


#TODO Make this a test
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