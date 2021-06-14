
####################################################
# Cell Division
####################################################

#This is so easy with a graph 
function CellDivision!(CPM::CellPotts, σ::Int)

    #Find all the nodes with the same ID (σ) as input
    cellIdx = findall(isequal(σ), CPM.graph.σ)

    #Use Metis to evenly partition the cell into two subcells
    newCellsIdx = Metis.partition(CPM.graph.network[cellIdx], 2) #Vector of ones and twos
    newCellNodeIdx = cellIdx[newCellsIdx .== 2]
    
    #Now we just need to update the model to account for the new cell

    #Update cell attributes
    #ID
        push!(CPM.cell.ids, CPM.cell.ids[end]+1)
    #Current volumes
        CPM.cell.volumes[σ] = sum(x->x==1, newCellsIdx)
        push!(CPM.cell.volumes, sum(x->x==2, newCellsIdx))

    #desired volumes
        push!(CPM.cell.desiredVolumes, CPM.cell.desiredVolumes[σ])
        
    #Types
        push!(CPM.cell.types, CPM.cell.types[σ])

    #Update graph attributes
    #σ
        CPM.graph.σ[newCellNodeIdx] .= CPM.cell.ids[end]
    #τ
        CPM.graph.τ[newCellNodeIdx] .= CPM.cell.types[end]
    #isArticulation
        UpdateConnections!(CPM.graph)

    #Penalties???   

    #Global CellPotts attributes
    #energy
        CPM.energy = sum([f(CPM) for f in CPM.M.penalties])
    #visual
        CPM.visual[newCellNodeIdx] .= CPM.cell.ids[end]

    return nothing
end