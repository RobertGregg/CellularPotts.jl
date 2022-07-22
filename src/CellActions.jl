#TODO either add removeCell() or condition for division to not select dead cells

####################################################
# Cell Division
####################################################

#This is so easy with a graph 
function CellDivision!(cpm::CellPotts, σ::Int)

    #Find all the nodes with the same ID (σ) as input
    cellNodeIDs = findall(isequal(σ), cpm.space.nodeIDs)

    #Use Metis to evenly partition the cell into two subcells
    nodePartition = Metis.partition(cpm.space[cellNodeIDs], 2) #Vector of ones and twos
    newCellNodeIDs = cellNodeIDs[nodePartition .== 2]
    
    #Now we just need to update the model to account for the new cell

    #Update old cell
        cpm.currentState.volumes[σ] = count(isequal(1), nodePartition)
    #Add new cell
        cell = cpm.currentState[σ]
        cell.cellIDs[1] = maximum(cpm.currentState.cellIDs) + 1 #now have n+1 cells
        cell.volumes[1] = count(isequal(2), nodePartition)
        addNewCell(cpm.currentState, cell)


    #Update graph attributes
        cpm.space.nodeIDs[newCellNodeIDs] .= cpm.currentState.cellIDs[end]

    #Global CellPotts attributes
    #visual
        cpm.visual[newCellNodeIDs] .= cpm.currentState.typeIDs[end]

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