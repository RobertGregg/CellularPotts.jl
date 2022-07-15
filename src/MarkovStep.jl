####################################################
# Metropolis–Hasting Step
####################################################

function MHStep!(cpm::CellPotts)

    #unpack current step structure update
    step = cpm.step


    #Loop through until a good source target is found
    searching = true
    while searching
        #Pick a random location on the graph
        step.sourceNode = rand(1:nv(cpm.space))
        #What cell does it belong to?
        step.sourceCellID = cpm.space.nodeIDs[step.sourceNode]

        #Get all of the unique cell IDs neighboring this Node
        step.neighborNodes = neighbors(cpm.space, step.sourceNode)

        #Choose a target
        step.targetCellID = cpm.space.nodeIDs[rand(step.neighborNodes)]

        #Some checks before attempting a flip
            #In the middle of a cell
            inMiddle = checkMiddle(cpm, step)
            #target is the same as source cell 
            isSource = step.targetCellID == step.sourceCellID

        #if all checks pass, attempt flip
        if !(inMiddle | isSource) 
            searching = false
        end
    end    



    #Calculate the change in energy when source node is modified
    #ΔH =  sum(f(cpm, step) for f in cpm.parameters.penalties)
    ΔH =  applyPenalties(cpm)

    #Calculate an acceptance ratio
    acceptRatio = min(1.0,exp(-ΔH/cpm.temperature))


    if rand() < acceptRatio #If the acceptance ratio is large enough
        #Need to update all cell and graph properties
        #---Cell properties---

        #Update cell volumes
        cpm.currentState.volumes[step.sourceCellID] -= 1
        cpm.currentState.volumes[step.targetCellID] += 1

        #TODO update cell perimeters

        #---Graph properties---

        #Cell IDs
        cpm.space.nodeIDs[step.sourceNode] = step.targetCellID
        cpm.space.nodeTypes[step.sourceNode] = cpm.currentState.names[step.targetCellID]
        
        #---Overall properties---
        #Update visual
        cpm.visual[step.sourceNode] = step.targetCellID
        cpm.step.stepCounter += 1

    end

    return nothing
end

function checkMiddle(cpm, step)
    
    for neighbor in step.neighborNodes
        if step.sourceCellID ≠ cpm.space.nodeIDs[neighbor]
            return false
        end
    end

    return true
end

function applyPenalties(cpm)
    ΔH = 0

    for penalty in cpm.penalties
        ΔH += penalty(cpm)
    end

    return ΔH
end