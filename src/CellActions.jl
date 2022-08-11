####################################################
# Cell Division
####################################################

#TODO User needs to specify what happens to custom properties
function CellDivision!(cpm::CellPotts, σ::Int)

    #Find all the nodes with the same ID (σ) as input
    cellNodeIDs = findall(isequal(σ), cpm.space.nodeIDs[:])

    #Use Metis to evenly partition the cell into two subcells
    nodePartition = Metis.partition(cpm.space[cellNodeIDs], 2) #Vector of ones and twos
    newCellNodeIDs = cellNodeIDs[nodePartition .== 2]
    
    #Now we just need to update the model to account for the new cell

    #Update old cell
    cpm.currentState.volumes[σ] = count(isequal(1), nodePartition)

    #Add new cell
    cellcopy = cpm.currentState[σ]
    addNewCell(cpm.currentState, cellcopy)

    #Update new cell attributes
    cpm.currentState.cellIDs[end] = maximum(cpm.currentState.cellIDs) + 1
    cpm.currentState.volumes[end] = count(isequal(2), nodePartition)


    #Update graph attributes
        cpm.space.nodeIDs[newCellNodeIDs] .= cpm.currentState.cellIDs[end]

    return nothing
end


####################################################
# Cell Death
####################################################


#Just set the desired volume to zero?
function CellDeath!(cpm::CellPotts, σ::Int)
    cpm.currentState.desiredVolumes[σ] = 0
    return nothing
end