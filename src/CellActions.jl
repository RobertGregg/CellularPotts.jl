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

    #Add new cell
    cellcopy = cpm.state[σ]
    addnewcell(cpm.state, cellcopy)

    #Newly created cell ID
    σnew = countcells(cpm)

    #Update old cell
    cpm.state.volumes[σ] = count(isequal(1), nodePartition)

    #Update new cell attributes
    cpm.state.cellIDs[σnew] = σnew
    cpm.state.volumes[σnew] = count(isequal(2), nodePartition)

    #Update space attributes
    cpm.space.nodeIDs[newCellNodeIDs] .= cpm.state.cellIDs[σnew]

    #Recalculate the perimeters as well
    cpm.state.perimeters[σ] = calcuatePerimeter(cpm, σ)
    cpm.state.perimeters[σnew] = calcuatePerimeter(cpm, σnew)



    #Also need to record all of the updates
    if cpm.record
        #Old Cell
        push!(cpm.history[:state][:volumes], [cpm.step.stepCounter, σ, cpm.state.volumes[σ]])
        push!(cpm.history[:state][:perimeters], [cpm.step.stepCounter, σ, cpm.state.perimeters[σ]])

        #New Cell
        for property in keys(cpm.state)
            push!(cpm.history[:state][property], [cpm.step.stepCounter, σnew, getproperty(cpm.state,property)[σnew]])
        end

        #Space
        for i in newCellNodeIDs
            push!(cpm.history[:space][:nodeIDs], [cpm.step.stepCounter, i, σnew])
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