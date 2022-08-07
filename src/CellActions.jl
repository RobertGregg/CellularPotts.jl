####################################################
# Cell Division
####################################################

#TODO Syntax now doesn't match b/s we're using Tables
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
        cell = cpm.currentState[σ]
        cell.cellIDs = maximum(cpm.currentState.cellIDs) + 1 #now have n+1 cells
        cell.volumes = count(isequal(2), nodePartition)
        addNewCell(cpm.currentState, cell)


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