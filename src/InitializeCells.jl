
#Presumably I will add more ways to initialize cells...

####################################################
# Functions to set initial cell locations
####################################################


#TODO Allow image input
#need two inputs:
    #1) image colored by cell ID
    #2) image colored by cell type


function positionCellsRandom!(cpm::CellPotts)

    #Generate random centers for each cell
    #TODO Use Poisson Disc sampling?
    centers = randperm(nv(cpm.space))[1:countcells(cpm)]

    cellMembership = growcells(cpm, centers)

    updateCellMembership!(cpm, cellMembership)
    
    return nothing
end


function positionCells!(cpm::CellPotts)

    #Convert positions to linear LinearIndices
    centers = [LinearIndices(size(cpm.space))[i...] for i in cpm.state.positions]

    cellMembership = growcells(cpm, centers)


    updateCellMembership!(cpm, cellMembership)

    return nothing
end


function growcells(cpm::CellPotts, centers)

    #Unpack some needed variables
    cellIDs = cpm.state.cellIDs
    volumes = cpm.state.desiredVolumes
    space = cpm.space

    numCells = countcells(cpm)

    #Using BFS in to grow each cell so need a queue
    queues = [[node] for node in centers]

    cellMembership = zeros(Int, size(space))
    cellMembership[centers] .= 1:numCells

    #Keep track of cells that are still growing as well as their current volume
    cellsStillGrowing = collect(1:numCells)
    cellSizes = ones(Int, numCells)

    while !isempty(cellsStillGrowing)

        #Pick a random growing cell
        i = rand(cellsStillGrowing)
        q = queues[i]

        #If the queue is empty, the cell has no room to grow and is removed 
        if isempty(q)
            filter!(!isequal(i), cellsStillGrowing)
            continue
        end

        node = popfirst!(q)
        for neighbor in neighbors(space, node)
            
            #Has that neighboring site been occupied?
            if cellMembership[neighbor] ≠ 0
                continue
            end

            #If not, add it to the cell and increase cell size
            cellMembership[neighbor] = cellIDs[i]
            cellSizes[i] += 1

            #If the cell is the correct size remove it from growing cells
            if cellSizes[i] == volumes[i]
                filter!(!isequal(i), cellsStillGrowing)
            end

            #Add to the cell queue
            push!(q, neighbor)

        end
    end

    return cellMembership
end

#Once the cells are positioned, update the model
function updateCellMembership!(cpm, cellMembership)
    
    #Update the network with the new cell locations
    for (i, cellID) in enumerate(cellMembership)
        if cellID ≠ 0
            cpm.space.nodeIDs[i] = cellID
            cpm.space.nodeTypes[i] = cpm.state.typeIDs[cellID]

            #Update the cell summary volumes
            cpm.state.volumes[cellID] += 1
        end
    end
    
    #Also update the cell perimeters
    #Can't be done in previous loop b/c not all nodeIDs are updated
    for (i, cellID) in enumerate(cellMembership)
        if cellID ≠ 0
            for n in neighbors(cpm.space, i)
                if cpm.space.nodeIDs[n] ≠ cellID
                    cpm.state.perimeters[cellID] += 1
                end
            end
        end
    end
end