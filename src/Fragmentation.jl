#Returns true is deleting the target node would disconnect the graph
#Uses a breadth first search to determine if all vertices can be reached after node removal
#Not faster than Tajan's algorithm but much simplier

#TODO Need to find the allocations
function isfragmented(cpm::CellPotts)

    #Unpack cpm and reset state (no allocations right?)
    g = cpm.space
    target = cpm.step.target
    F = cpm.fragment

    #TODO Not checking Medium mean cells can enclose bits of medium
    #No need to check Medium for articulation points  
    if iszero(target.id)
        return false
    end
    
    resetArticulation!(F)
    
    #find a neighbor to the target node that has the same ID as the target
    startID = findfirst(isequal(target.id), g.nodeIDs[i] for i in target.neighbors)
    startNode = target.neighbors[startID]

    #Mark the target and neighbor node as visited
    F.visited[target.node] = true
    F.visited[startNode] = true

    #Add the target neighbor to the queue
    push!(F.queue, startNode)

    #Breath first search through the cell
    while !isempty(F.queue)
        currentNode = popfirst!(F.queue)

        #Loop through neighbors in target cell and not equal target ID
        for neighbor in cellneighbors(g,currentNode)
            if !F.visited[neighbor]
                F.visited[neighbor] = true
                push!(F.queue, neighbor)
            end

            # #TODO Need to check this, but if we reach all the neighbors of target node then I think we can stop
            if all(view(F.visited, target.neighbors))
                return false
            end
        end
    end

    #TODO Use state volume? This won't work for Medium.
    #If we visited all the nodes (except target node) then cell is still connected
    if count(isequal(target.id), g.nodeIDs) == count(F.visited)
        return false
    end

    return true
end



#Special case for 2D, planar spaces
#Count the number of vertices, edges, and faces changed when removing "node" from the graph
#V-E+F == 2 for all planar connected graphs, therefore we should see no change when removing "node"
function isfragmented(cpm::CellPotts{T,4,2,S,U}) where {T,S,U}
    
    g = cpm.space
    node = cpm.step.target.node
    id = cpm.step.target.id

    V = 1
    
    E = count(g.nodeIDs[n] == id for n in neighbors(g,node))
    
    F = 0
    for (x,y) in allpairs(cellneighbors(g,node))
        if count(n -> n ∈ cellneighbors(g,x), cellneighbors(g,y)) == 2
            F += 1
        end
    end

    return iszero(V-E+F) ? false : true
end
    
#TODO 3D version is possible
# V-E+F-B=1 where B is the number of "boxes"