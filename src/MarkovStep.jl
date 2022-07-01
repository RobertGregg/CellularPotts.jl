####################################################
# Metropolis–Hasting Step
####################################################

function MHStep!(cpm::CellPotts)

    #unpack current step structure update
    stepInfo = cpm.stepInfo


    #Loop through until a good source target is found

    searching = true
    while searching
        #Pick a random location on the graph
        stepInfo.sourceNode = rand(1:nv(cpm.graph))
        #What cell does it belong to?
        stepInfo.sourceCell = cpm.graph.nodeIDs[stepInfo.sourceNode]

        #Get all of the unique cell IDs neighboring this Node
        stepInfo.sourceNodeNeighbors = neighbors(cpm.graph, stepInfo.sourceNode)
        stepInfo.possibleCellTargets = cpm.graph.nodeIDs[stepInfo.sourceNodeNeighbors]
        unique!(stepInfo.possibleCellTargets)
        #Choose a target
        stepInfo.targetCell = rand(stepInfo.possibleCellTargets)

        #Collect the target and source into a vector
        stepInfo.sourceTargetCell[1] = stepInfo.sourceCell
        stepInfo.sourceTargetCell[2] = stepInfo.targetCell

        #Some checks before attempting a flip
            #In the middle of a cell
            inMiddle = length(stepInfo.possibleCellTargets) == 1
            #target is the same as source cell 
            isSource = stepInfo.targetCell == stepInfo.sourceCell

        #if all checks pass, attempt flip
        if !(inMiddle | isSource) 
            searching = false
        end
    end    



    #Calculate the change in energy when source node is modified
    ΔH =  sum(f(cpm, stepInfo) for f in cpm.parameters.penalties)

    #Calculate an acceptance ratio
    acceptRatio = min(1.0,exp(-ΔH/cpm.parameters.temperature))



    if rand() < acceptRatio #If the acceptance ratio is large enough

        #Need to update all cell and graph properties
        #---Cell properties---

        #Update cell volumes
        cpm.cells.volumes[stepInfo.sourceCell] -= 1
        cpm.cells.volumes[stepInfo.targetCell] += 1


        #---Graph properties---

        #Cell IDs
        cpm.graph.nodeIDs[stepInfo.sourceNode] = stepInfo.targetCell
        cpm.graph.nodeTypes[stepInfo.sourceNode] = cpm.cells.names[stepInfo.targetCell]
        
        #---Overall properties---
        #Update energy
        cpm.energy += ΔH
        #Update visual
        cpm.visual[stepInfo.sourceNode] = stepInfo.targetCell
        cpm.stepCounter += 1

    end

    return nothing
end