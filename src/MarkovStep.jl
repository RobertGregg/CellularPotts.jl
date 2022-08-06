####################################################
# Helper Functions
####################################################

function calculateΔH(cpm::CellPotts)
    ΔH = 0

    #loop though indices to exploit known union types
    for i in eachindex(cpm.penalties)
        ΔH += addPenalty!(cpm, cpm.penalties[i])
    end

    return ΔH
end

####################################################
# Metropolis–Hasting Step
####################################################

function MHStep!(cpm::CellPotts)
    
    #Pick a random location on the graph
    cpm.step.sourceNode = rand(1:nv(cpm.space))
    #What cell does it belong to?
    cpm.step.sourceCellID = cpm.space.nodeIDs[cpm.step.sourceNode]

    #Get all cell IDs neighboring this Node
    cpm.step.sourceNeighborNodes = neighbors(cpm.space, cpm.step.sourceNode)

    #Choose a target
    cpm.step.targetNode = rand(cpm.step.sourceNeighborNodes)
    cpm.step.targetCellID = cpm.space.nodeIDs[cpm.step.targetNode]
    
    #Some checks before attempting a flip

    #target and source same cell? 
    if cpm.step.targetCellID == cpm.step.sourceCellID
        return nothing
    end

    #sourceCell surrounded by only sourceCells?
    if all(isequal(cpm.step.sourceCellID, cpm.space.nodeIDs[n]) for n in cpm.step.sourceNeighborNodes)
        return nothing
    end

    #Test if the target cell is an articulation point
    if cpm.step.targetNode ∈ cpm.getArticulation(cpm.space, cpm.step.targetCellID)
        return nothing
    end

    #Get all cell nodes neighboring for target node
    cpm.step.targetNeighborNodes = neighbors(cpm.space, cpm.step.targetNode)

    #Calculate the change in energy when target node is modified
    ΔH =  calculateΔH(cpm)


    #Calculate an acceptance ratio
    acceptRatio = min(1.0,exp(-ΔH/cpm.temperature))


    if rand() < acceptRatio #If the acceptance ratio is large enough
        #Need to update all cell and graph properties
        #---Graph properties---

        #Cell IDs
        cpm.space.nodeIDs[cpm.step.targetNode] = cpm.step.sourceCellID
        cpm.space.nodeTypes[cpm.step.targetNode] = cpm.currentState.typeIDs[cpm.step.sourceCellID]

        #---Cell properties---
        for i in eachindex(cpm.penalties)
            updateMHStep!(cpm, cpm.penalties[i])
        end
        
    end

    return nothing
end


####################################################
# Model Step
####################################################

function ModelStep!(cpm::CellPotts)

    #Repeat MHStep! for each pixel in the model
    for _ in 1:prod(cpm.space.gridSize)
        MHStep!(cpm)
    end

    #Increment the step counter
    cpm.step.stepCounter += 1

    #Model updates after MCS
    for i in eachindex(cpm.penalties)
        updateModelStep!(cpm, cpm.penalties[i])
    end

    return nothing
end


####################################################
# How penalties update after each MHStep!
####################################################

#Default to nothing
updateMHStep!(cpm::CellPotts, penalty::Penalty) = nothing

function updateMHStep!(cpm::CellPotts, VP::VolumePenalty)
    #Update cell volumes
    cpm.currentState.volumes[cpm.step.sourceCellID] += 1
    cpm.currentState.volumes[cpm.step.targetCellID] -= 1
    return nothing
end

function updateMHStep!(cpm::CellPotts, PP::PerimeterPenalty)
    #Update cell perimeters
    cpm.currentState.perimeters[cpm.step.sourceCellID] += PP.Δpᵢ
    cpm.currentState.perimeters[cpm.step.targetCellID] -= PP.Δpⱼ
    return nothing
end


function updateMHStep!(cpm::CellPotts, MP::MigrationPenalty)

    τ = cpm.currentState.typeIDs[cpm.step.targetCellID]
    
    #Do not update cells with λ==0
    if iszero(MP.λ[τ])
        MP.nodeMemory[cpm.step.targetNode] = 0
    else
        MP.nodeMemory[cpm.step.targetNode] = MP.maxAct
    end
    
    return nothing
end

####################################################
# How penalties update after each ModelStep!
####################################################

#Default to nothing
updateModelStep!(cpm::CellPotts, penalty::Penalty) = nothing

function updateModelStep!(cpm::CellPotts, MP::MigrationPenalty)
    
    #Remove zeros from Medium copying into cells
    dropzeros!(MP.nodeMemory)

    #Reduce the node Memory by 1 and remove zeros
    MP.nodeMemory.nzval .-= 1
    dropzeros!(MP.nodeMemory)
    
    return nothing
end