####################################################
# Metropolis–Hasting Step
####################################################

function MHStep!(cpm::CellPotts)

    #unpack current step structure
    step = cpm.step
    

    #Pick a random location on the graph
    step.sourceNode = rand(1:nv(cpm.space))
    #What cell does it belong to?
    step.sourceCellID = cpm.space.nodeIDs[step.sourceNode]

    #Get all of the unique cell IDs neighboring this Node
    step.neighborNodes = neighbors(cpm.space, step.sourceNode)

    #Choose a target
    step.targetNode = rand(step.neighborNodes)
    step.targetCellID = cpm.space.nodeIDs[step.targetNode]

    #Some checks before attempting a flip
    #In the middle of a cell
    if all(isequal(step.sourceCellID, cpm.space.nodeIDs[n]) for n in step.neighborNodes)
        return nothing
    end

    #target is the same as source cell 
    if step.targetCellID == step.sourceCellID
        return nothing
    end



    #Calculate the change in energy when target node is modified
    #ΔH =  sum(f(cpm, step) for f in cpm.parameters.penalties)
    ΔH =  applyPenalties(cpm)

    #Calculate an acceptance ratio
    acceptRatio = min(1.0,exp(-ΔH/cpm.temperature))


    if rand() < acceptRatio #If the acceptance ratio is large enough
        #Need to update all cell and graph properties
        #---Cell properties---

        for i in eachindex(cpm.penalties)
            updateStep!(cpm, step, cpm.penalties[i])
        end

        #---Graph properties---

        #Cell IDs
        cpm.space.nodeIDs[step.sourceNode] = step.targetCellID
        cpm.space.nodeTypes[step.sourceNode] = cpm.currentState.names[step.targetCellID]
        
        #TODO Add articulation point updater 
        
        #---Overall properties---
        #Update visual
        cpm.visual[step.sourceNode] = cpm.currentState.typeIDs[step.targetCellID]
    end

    return nothing
end


####################################################
# Model Step
####################################################
function ModelStep!(cpm::CellPotts)
    for _ in 1:prod(cpm.space.gridSize)
        MHStep!(cpm)
    end
    cpm.step.stepCounter += 1

    return nothing
end


updateStep!(cpm::CellPotts, step::MHStepInfo, AP::AdhesionPenalty) = nothing

function updateStep!(cpm::CellPotts, step::MHStepInfo, VP::VolumePenalty)
    #Update cell volumes
    cpm.currentState.volumes[step.sourceCellID] -= 1
    cpm.currentState.volumes[step.targetCellID] += 1
    return nothing
end

function updateStep!(cpm::CellPotts, step::MHStepInfo, PP::PerimeterPenalty)
    #Update cell perimeters
    cpm.currentState.perimeters[step.sourceCellID] -= PP.Δpᵢ
    cpm.currentState.perimeters[step.targetCellID] += PP.Δpⱼ
    return nothing
end


function updateStep!(cpm::CellPotts, step::MHStepInfo, MP::MigrationPenalty)
    
    #Reduce the node Mo=emory by 1 and remove zeros
    MP.nodeMemory.nzval .-= 1
    dropzeros!(MP.nodeMemory)
    
    #Update cell volumes
    MP.nodeMemory[cpm.step.targetNode] += MP.maxAct
    return nothing
end

function applyPenalties(cpm)
    ΔH = 0

    for i in eachindex(cpm.penalties)
        ΔH += addPenalty!(cpm, cpm.penalties[i])
    end

    return ΔH
end