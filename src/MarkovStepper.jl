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

    source = cpm.step.source
    target = cpm.step.target

    #Pick a random location on the graph
    source.node = rand(1:nv(cpm.space))
    #What cell does it belong to?
    source.id = cpm.space.nodeIDs[source.node]
    source.type = cpm.space.nodeTypes[source.node]

    #Get all cell IDs neighboring this Node
    source.neighbors = neighbors(cpm.space, source.node)

    #Some space locations are forbidden (see TightSpaces example) 
    if isempty(source.neighbors)
        return nothing
    end

    #Choose a target
    target.node = rand(source.neighbors)
    target.id = cpm.space.nodeIDs[target.node]
    target.type = cpm.space.nodeTypes[target.node]
    
    #Some checks before attempting a flip

    #target and source same cell? 
    if target.id == source.id
        return nothing
    end

    #sourceCell surrounded by only sourceCells?
    if all(isequal(source.id, cpm.space.nodeIDs[n]) for n in source.neighbors)
        return nothing
    end

    
    #Get all cell nodes neighboring target node
    target.neighbors = neighbors(cpm.space, target.node)
    
    #Calculate the change in energy when target node is modified
    ΔH =  calculateΔH(cpm)
    
    
    #Calculate an acceptance ratio
    acceptRatio = min(1.0,exp(-ΔH/cpm.temperature))
    
    
    if rand() < acceptRatio #If the acceptance ratio is large enough

        #Moved into accept loop b/c computationally intensive
        if cpm.checkArticulation && isfragmented(cpm)
            return nothing
        end

        if cpm.recordHistory
            updateHistory!(cpm)
        end

        #Cell IDs
        cpm.space.nodeIDs[target.node] = source.id
        cpm.space.nodeTypes[target.node] = source.type


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
    cpm.step.counter += 1

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
    cpm.state.volumes[cpm.step.source.id] += 1
    cpm.state.volumes[cpm.step.target.id] -= 1

    return nothing
end

function updateMHStep!(cpm::CellPotts, PP::PerimeterPenalty)
    #Update cell perimeters
    cpm.state.perimeters[cpm.step.source.id] += PP.Δpᵢ
    cpm.state.perimeters[cpm.step.target.id] -= PP.Δpⱼ

    return nothing
end


function updateMHStep!(cpm::CellPotts, MP::MigrationPenalty)

    τ = cpm.state.typeIDs[cpm.step.source.id]
    
    #Do not update cells with λ==0
    if iszero(MP.λ[τ])
        MP.nodeMemory[cpm.step.target.node] = 0
    else
        MP.nodeMemory[cpm.step.target.node] = MP.maxAct
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