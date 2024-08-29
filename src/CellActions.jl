####################################################
# Cell Division
####################################################

#TODO User needs to specify what happens to custom properties
function CellDivision!(cpm::CellPotts, oldID::Int)

    #Find all the nodes with the same ID (σ) as input
    cellNodeIDs = findall(isequal(oldID), cpm.space.nodeIDs[:])

    #Use Metis to evenly partition the cell into two subcells
    nodePartition = Metis.partition(cpm.space[cellNodeIDs], 2) #Vector of ones and twos
    newCellNodeIDs = cellNodeIDs[nodePartition .== 2]
    

    #Now we just need to update the model to account for the new cell

    #Add new cell (copies properties from parent cell)
    #Updates cellID automatically
    push!(cpm.state, cpm.state[oldID])

    #Newly created cell ID
    newID = cpm.state.cellIDs[end]

    #Update old cell
    cpm.state.volumes[oldID] = count(isequal(1), nodePartition)

    #Update new cell attributes
    cpm.state.volumes[newID] = count(isequal(2), nodePartition)

    #Update space attributes
    cpm.space.nodeIDs[newCellNodeIDs] .= cpm.state.cellIDs[newID]

    #Recalculate the perimeters as well
    cpm.state.perimeters[oldID] = calcuatePerimeter(cpm, oldID)
    cpm.state.perimeters[newID] = calcuatePerimeter(cpm, newID)



    #Also need to record all of the updates
    if cpm.recordHistory
        for i in newCellNodeIDs
            updateHistory!(cpm, cpm.step.counter, i, newID, cpm.state.cellIDs[newID])
        end
    end

    return nothing
end


####################################################
# Cell Death
####################################################


#Just set the desired volume to zero?
function CellDeath!(cpm::CellPotts, σ::Int)
    cpm.state.desiredVolumes[σ] = 0
    return nothing
end