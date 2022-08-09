
#Presumably I will add more ways to initialize cells...

####################################################
# Functions to set initial cell locations
####################################################


#TODO Allow image input
#need two inputs:
    #1) image colored by cell ID
    #2) image colored by cell type


#Once the cells are positioned, update the model
function updateCellMembership!(cpm, cellMembership)
    
    #Update the network with the new cell locations
    for (i, cellID) in enumerate(cellMembership)
        if cellID ≠ 0
            cpm.space.nodeIDs[i] = cellID
            cpm.space.nodeTypes[i] = cpm.currentState.typeIDs[cellID]

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
end


function positionCellsRandom!(cpm::CellPotts{N,T,V}) where {N,T,V}

    #Unpack the cell space
    space = cpm.space

    #initialize matrix of cell IDs (σ)
    cellMembership = zeros(T, space.gridSize)
    
    #Find the center node of the entire graph
    centerIdx = CartesianIndex(size(space).÷2)
    nodeIdx = LinearIndices(size(space))
    centerNode = nodeIdx[centerIdx]

    #Determine how far nodes are from the center
    nodeDis = gdistances(space, centerNode)

    #Get a sorted permutation of the node distances
    sortedDis = sortperm(nodeDis)

    #How many nodes need to be initialized?
    totalNodes = sum(cpm.currentState.desiredVolumes)

    #Assign 1:totalNodes to be filled with cells and the rest medium
    networkIdx = sortedDis[1:totalNodes] 

    #Partition the identified nodes by the number of cells needed
    if countcells(cpm) == 1 #There is only one cell (no need to partition)
        cellMembership[networkIdx] .= 1
    else
        cellMembership[networkIdx] = Metis.partition(cpm.space[networkIdx], countcells(cpm))
    end

    #Shuffle around the cells IDs so similar cells types are not clumped together
    replace!(cellMembership, (1:countcells(cpm) .=> shuffle(1:countcells(cpm)))...)

    updateCellMembership!(cpm, cellMembership)
    
    return nothing
end


function positionCells!(cpm::CellPotts{N,T,V}) where {N,T,V}

    #Unpack the initial state and parameters
    space = cpm.space
    currState = cpm.currentState

    #initialize matrix of cell IDs (σ)
    cellMembership = zeros(T, space.gridSize)

    for cell in currState

        if iszero(cell.cellIDs)
            continue
        end

        #Determine how far nodes are from the position
        node = LinearIndices(size(space))[cell.positions...]
        
        politeBFS(cellMembership, space, cell.cellIDs, cell.desiredVolumes, node)
    end


    updateCellMembership!(cpm, cellMembership)

    return nothing
end


function politeBFS(cellMembership, space, id, vol, node)

    queue = [node]
    cellMembership[node] = id
    cellSize = 1

    while !isempty(queue)
        nextNode = popfirst!(queue)
        for neighbor in neighbors(space,nextNode)
            if cellMembership[neighbor] == 0
                cellMembership[neighbor] = id
                cellSize += 1
                if cellSize == vol
                    return nothing
                end
                push!(queue, neighbor)
            end
        end
    end

    return nothing
end
