
#Presumably I will add more ways to initialize cells...

####################################################
# Functions to set initial cell locations
####################################################


#TODO Allow image input
#TODO Allow position input


function positionCellsRandom!(cpm::CellPotts{N,T,V}) where {N,T,V}

    #Unpack the initial state and parameters
    space = cpm.space

    #initialize matrix of cell IDs (σ)
    cellMembership = zeros(T, space.gridSize)
    
    #Find the center node of the entire graph
    centerIdx = CartesianIndex(space.gridSize.÷2)
    nodeIdx = LinearIndices(space.gridSize)
    centerNode = nodeIdx[centerIdx]

    #Determine how far nodes are from the center
    nodeDis = gdistances(space, centerNode)

    #How many nodes need to be initialized?
    totalNodes = sum(cpm.currentState.desiredVolumes)

    #Get a sorted permutation of the distance
    sortedDis = sortperm(nodeDis)

    #Assign 1:totalNodes to be filled with cells and the rest medium
    networkIdx = sortedDis[1:totalNodes] 

    #Partition the identified nodes by the number of cells needed
    if countCellTypes(cpm) == 1 #There is only one cell (no need to partition)
        cellMembership[networkIdx] .= 1
    else
        cellMembership[networkIdx] = Metis.partition(cpm.space[networkIdx], countCells(cpm))
    end

    #Update the network with the new cell locations
    for (i, cellID) in enumerate(cellMembership)
        if cellID ≠ 0
            cpm.space.nodeIDs[i] = cellID
            cpm.space.nodeTypes[i] = cpm.currentState.names[cellID]

            #Update the cell summary volumes
            cpm.currentState.volumes[cellID] += 1
        end
    end
    
    #Also update the cell perimeters
    #Can't be done in previous loop b/c not all nodeIDs are updated
    for (i, cellID) in enumerate(cellMembership)
        if cellID ≠ 0
            for n in neighbors(cpm.space, i)
                if cpm.space.nodeIDs[n] ≠ cellID
                    cpm.currentState.perimeters[cellID] += 1
                end
            end
        end
    end

    #Fill the the array for the visual
    cpm.visual = cellMembership

    return nothing
end


function positionCells!(cpm::CellPotts{N,T,V}) where {N,T,V}

    #Unpack the initial state and parameters
    space = cpm.space
    currState = cpm.currentState

    #initialize matrix of cell IDs (σ)
    cellMembership = zeros(T, space.gridSize)

    for (id, vol, pos) in zip(currState.cellIDs, currState.desiredVolumes, currState.positions)

        #Determine how far nodes are from the position
        node = LinearIndices(space.gridSize)[pos...]
        nodeDis = gdistances(space, node)

        #Get a sorted permutation of the distance
        sortedDis = sortperm(nodeDis)

        #Assign 1:totalNodes to be filled with cells and the rest medium
        networkIdx = sortedDis[1:vol]

        cellMembership[networkIdx] .= id
    end


    #Update the network with the new cell locations
    for (i, cellID) in enumerate(cellMembership)
        if cellID ≠ 0
            cpm.space.nodeIDs[i] = cellID
            cpm.space.nodeTypes[i] = cpm.currentState.names[cellID]

            #Update the cell summary volumes
            cpm.currentState.volumes[cellID] += 1
        end
    end
    
    #Also update the cell perimeters
    #Can't be done in previous loop b/c not all nodeIDs are updated
    for (i, cellID) in enumerate(cellMembership)
        if cellID ≠ 0
            for n in neighbors(space, i)
                if space.nodeIDs[n] ≠ cellID
                    cpm.currentState.perimeters[cellID] += 1
                end
            end
        end
    end

    #Fill the the array for the visual
    cpm.visual = cellMembership

    return nothing
end

