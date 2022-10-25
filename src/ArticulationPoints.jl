####################################################
# Utility to get articulation points
####################################################

#=
NOTE: Graphs.jl has an articulation() function which is very good. I went through the trouble of creating this utility to avoid allocations because articulation checks are done after each MHStep. Depending on the space size this gives a 1.2x to 5x speedup.
=#

"""
    ArticulationUtility(n::Int)
Create a new `ArticulationUtility` to hold variables used in articulation point calculation.

This data structure contains the following fields:
 - discoveryTime`::Vector{Int}`: The time step where the node was discovered
 - low`::Vector{Int}`: Minimum discovery time when you account for back edges
 - parent`::Vector{Int}`: Parent node (i.e. how we got to this node)
 - visited`::BitVector`: Track if node has been visited
 - articulationPoints`::Vector{Int}`: List of articulation vertices
 - nodeIDs`::Vector{Int}`: Tracks the node IDs belonging to the current cell
 - time`::Int`: Increments by 1 every time a new node is discovered
 - cellID``::Int`: Cell ID (i.e. subgraph) we are testing for

Here `n`` is the number of vertices in the graph being tested
"""
mutable struct ArticulationUtility
    discoveryTime::Vector{Int}       #The time step where the node was discovered
    low::Vector{Int}                 #Minimum discovery time when you account for back edges   
    parent::Vector{Int}              #Parent node (i.e. how we got to this node)
    visited::BitVector               #Track if node has been visited
    articulationPoints::Vector{Int}  #List of articulation vertices
    nodeIDs::Vector{Int}             #Tracks the node IDs belonging to the current cell
    time::Int                        #Increments by 1 every time a new node is discovered
    cellID::Int                      #Cell ID (i.e. subgraph) we are testing for

    # Here n is the number of vertices in the graph
    function ArticulationUtility(n::Int)
        discoveryTime = zeros(Int,n)
        low = zeros(Int,n)
        parent = zeros(Int,n)

        visited = falses(n)
        articulationPoints = Int64[]
        sizehint!(articulationPoints, n) #max number of vertices

        nodeIDs = Int64[]
        sizehint!(nodeIDs, n) #max number of vertices

        time = 1
        cellID = 0

        return new(discoveryTime, low, parent, visited, articulationPoints, nodeIDs, time, cellID)
    end
end

#Set everything back to default values
function resetArticulationUtility(A::ArticulationUtility)

    A.discoveryTime[A.nodeIDs] .= 0
    A.low[A.nodeIDs] .= 0
    A.parent[A.nodeIDs] .= 0
    A.visited .= false
    empty!(A.nodeIDs)
    empty!(A.articulationPoints)
    A.time = 1
    A.cellID = 0

    return nothing
end

#Check if node i is a root node
isroot(A::ArticulationUtility, i::Int) = iszero(A.parent[i])

#Add node i to the list of articulation points
function pushcell!(A::ArticulationUtility, i::Int)
    if i ∉ A.articulationPoints
        push!(A.articulationPoints, i)
    end

    return nothing
end

#The ArticulationUtility object takes a space object and cellID as input
function (A::ArticulationUtility)(g::CellSpace, cellID::Int)

    resetArticulationUtility(A)
    
    #No need to check Medium for articulation points  
    if iszero(cellID)
        return A.articulationPoints
    end

    A.cellID = cellID

    #Loop through all nodes with current cellID
    @inbounds for i in Iterators.filter(x->g.nodeIDs[x] == cellID, vertices(g))
        push!(A.nodeIDs, i)
        if !A.visited[i]
            cellDFS!(A, g, i)
        end
    end

    return sort!(A.articulationPoints) #maybe use some kind of sorted set data structure?
end


#Recursively call depth first search to traverse graph
function cellDFS!(A::ArticulationUtility, g::CellSpace, i::Int)

    #Count children for current node
    children = 0

    #Mark the current node as visited
    A.visited[i] = true

    #Initialize discovery time and low value
    A.discoveryTime[i] = A.time
    A.low[i] = A.time

    #Increment time because we found a new node
    A.time += 1

    #Loop through adjacenct nodes
    #Ignore neighboring nodes from other cells
    @inbounds for v in Iterators.filter(x->g.nodeIDs[x] == A.cellID, neighbors(g,i))

        # If v is not visited yet, then make it a child of i in DFS tree
        if !A.visited[v]
            A.parent[v] = i
            children += 1
            #Recursively call DFS for this vertex
            cellDFS!(A, g, v)

            # Check if the subtree rooted with v has a connection to one of the ancestors of i
            A.low[i] = min(A.low[i], A.low[v])

            # i is an articulation point in following cases:
            # (1) i is root of DFS tree and has two or more children.
            if isroot(A,i) && children > 1
                pushcell!(A,i)
            end

            #(2) If i is not root and low value of one of its child is more than discovery value of i.
            if !isroot(A,i) && A.low[v] ≥ A.discoveryTime[i]
                pushcell!(A,i)
            end
        
        #Update low value of i for parent function calls
        elseif v ≠ A.parent[i]
            A.low[i] = min(A.low[i], A.discoveryTime[v])
        end
    end
end