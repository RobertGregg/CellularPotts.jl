
#Presumably I will add more ways to initialize cells...

####################################################
# Functions to set initial cell locations
####################################################

#Creates a cluster of cells in the middle of the grid
function GrowCells(gr::SimpleGraph, M::ModelParameters)

    #initialize matrix of cell indicies (σ)
    cellMembership = zeros(Int, M.graphDimension...)
    
    #Find the center node
    centerIdx = CartesianIndex(M.graphDimension.÷2)
    nodeIdx = LinearIndices(M.graphDimension)
    centerNode = nodeIdx[centerIdx]

    #Determine how far nodes are from the center
    nodeDis = gdistances(gr, centerNode)

    #How many nodes need to be initialized?
    totalNodes = M.cellCounts ⋅ M.cellVolumes #dot product to sum and multiply

    #Get a sorted permutation of the distance
    sortedDis = sortperm(nodeDis)

    #Assign 1:totalNodes to be filled with cells and the rest medium
    networkIdx = sortedDis[1:totalNodes] 

    #Partition the identified nodes by the number of cells needed
    if sum(M.cellCounts) == 1 #There is only one cell (no need to partition)
        cellMembership[networkIdx] .= 1
    else
        cellMembership[networkIdx] = Metis.partition(gr[networkIdx], sum(M.cellCounts))
    end

    return cellMembership
end
