
####################################################
# Cell Division
####################################################

#This is so easy with a graph 
function CellDivision!(CPM::CellPotts, σ::Int)

    cellIdx = findall(isequal(σ), CPM.graph.σ)

    #Vector of ones and twos
    newCellsIdx = Metis.partition(CPM.graph.network[cellIdx], 2)
    newCellNodeIdx = cellIdx[newCellsIdx .== 2]
    

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