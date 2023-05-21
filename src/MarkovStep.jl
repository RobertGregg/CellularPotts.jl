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

    #Reset the success flag
    cpm.step.success = false
    
    #Pick a random location on the graph
    cpm.step.sourceNode = rand(1:nv(cpm.space))
    #What cell does it belong to?
    cpm.step.sourceCellID = cpm.space.nodeIDs[cpm.step.sourceNode]

    #Get all cell IDs neighboring this Node
    cpm.step.sourceNeighborNodes = neighbors(cpm.space, cpm.step.sourceNode)

    #Some space locations are forbidden (see TightSpaces example) 
    if isempty(cpm.step.sourceNeighborNodes)
        return nothing
    end

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

    #TODO Add an option to toggle this check
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
        #Need to update properties

        #Cell IDs
        cpm.space.nodeIDs[cpm.step.targetNode] = cpm.step.sourceCellID
        cpm.space.nodeTypes[cpm.step.targetNode] = cpm.state.typeIDs[cpm.step.sourceCellID]

        if cpm.record
            updateHist!(cpm)
        end

        #---Cell properties---
        for i in eachindex(cpm.penalties)
            updateMHStep!(cpm, cpm.penalties[i])
        end

        #Finally update the success flag 
        cpm.step.success = true
    end

    return nothing
end


####################################################
# Model Step
####################################################

function ModelStep!(cpm::CellPotts)

    #Repeat MHStep! for each pixel in the model
    for _ in 1:nv(cpm.space)
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
    cpm.state.volumes[cpm.step.sourceCellID] += 1
    cpm.state.volumes[cpm.step.targetCellID] -= 1

    return nothing
end

function updateMHStep!(cpm::CellPotts, PP::PerimeterPenalty)
    #Update cell perimeters
    cpm.state.perimeters[cpm.step.sourceCellID] += PP.Δpᵢ
    cpm.state.perimeters[cpm.step.targetCellID] -= PP.Δpⱼ

    return nothing
end


function updateMHStep!(cpm::CellPotts, MP::MigrationPenalty)

    τ = cpm.state.typeIDs[cpm.step.sourceCellID]
    
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
    
    #Remove zeros from λ==0 cells copying into other cells
    dropzeros!(MP.nodeMemory)

    #Reduce the node memory by 1 and remove explicit zeros
    MP.nodeMemory.nzval .-= 1
    dropzeros!(MP.nodeMemory)
    
    return nothing
end