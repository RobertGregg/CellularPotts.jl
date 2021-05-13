####################################################
# Metropolis–Hasting Step
####################################################

function MHStep!(CPM::CellPotts)

    #Create a structure to hold the source and target information
    stepInfo = MHStepInfo() #should add as input, might save a little time

    #Loop through until a good source target is found
    searching = true
    while searching
        #Pick a random location on the graph
        stepInfo.sourceNode = rand(1:nv(CPM.graph.network))
        #What cell does it belong to?
        stepInfo.sourceCell = CPM.graph.σ[stepInfo.sourceNode]

        #Get all of the unique cell IDs neighboring this Node
        stepInfo.sourceNodeNeighbors = neighbors(CPM.graph.network, stepInfo.sourceNode)
        possibleCellTargets = unique(CPM.graph.σ[stepInfo.sourceNodeNeighbors])
        #Choose a target
        stepInfo.targetCell = rand(possibleCellTargets)

        #Some checks before attempting a flip
            #In the middle of a cell
            inMiddle = all(possibleCellTargets .== stepInfo.sourceCell)
            #will fragment the cell
            isArticulation = CPM.graph.isArticulation[stepInfo.sourceNode]
            #target is the same as source cell 
            isSource = stepInfo.targetCell == stepInfo.sourceCell

        #if all checks pass, attempt flip
        if !(inMiddle | isArticulation | isSource) 
            searching = false
        else
            CPM.stepCounter += 1
        end
    end    

    #Calculate the change in energy when source node is modified
    ΔH =  sum([f(CPM, stepInfo) for f in CPM.M.penalties])

    #Calculate an acceptance ratio
    acceptRatio = min(1.0,exp(-ΔH/CPM.M.temperature))

    if rand() < acceptRatio #If the acceptance ratio is large enough

        #Need to update all cell and graph properties
        #---Cell properties---

        #Update cell volumes
        CPM.cell.volumes[stepInfo.sourceCell] -= 1
        CPM.cell.volumes[stepInfo.targetCell] += 1


        #---Graph properties---

        #Cell IDs
        CPM.graph.σ[stepInfo.sourceNode] = stepInfo.targetCell
        CPM.graph.τ[stepInfo.sourceNode] = CPM.cell.types[stepInfo.targetCell]

        #Articulation points
        UpdateConnections!(CPM.graph)
        
        #---Overall properties---
        #Update energy
        CPM.energy += ΔH
        #Update visual
        CPM.visual[stepInfo.sourceNode] = stepInfo.targetCell

    end
    return nothing
end
