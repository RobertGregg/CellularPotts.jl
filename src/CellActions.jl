
####################################################
# Cell Division
####################################################

#This is so easy with a graph 
function CellDivision!(cpm::CellPotts, σ::Int)

    #Find all the nodes with the same ID (σ) as input
    cellIdx = findall(isequal(σ), cpm.graph.nodeIDs)

    #Use Metis to evenly partition the cell into two subcells
    newCellsIdx = Metis.partition(cpm.graph[cellIdx], 2) #Vector of ones and twos
    newCellNodeIdx = cellIdx[newCellsIdx .== 2]
    
    #Now we just need to update the model to account for the new cell

    #Update cell summary
        push!(cpm.cells.names, cpm.cells.names[σ]) 
        push!(cpm.cells.ids, cpm.cells.ids[end]+1) 
        #Current volumes
        cpm.cells.volumes[σ] = sum(x->x==1, newCellsIdx)
        push!(cpm.cells.volumes, sum(x->x==2, newCellsIdx))
        #Desired volumes 
        push!(cpm.cells.desiredVolumes, cpm.cells.desiredVolumes[σ])

    #Update graph attributes
        cpm.graph.nodeIDs[newCellNodeIdx] .= cpm.cells.ids[end]

    #Global CellPotts attributes
    #energy
        for penalty in cpm.parameters.penalties
            penalty(cpm)
        end
    #visual
        cpm.visual[newCellNodeIdx] .= cpm.cells.ids[end]

    return nothing
end


####################################################
# Cell Death
####################################################


#Just set the desired volume to zero?
function CellDeath!(cpm::CellPotts, σ::Int)
    cpm.cells.desiredVolumes[σ] = 0
    return nothing
end