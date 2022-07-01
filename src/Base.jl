include("GraphStructure.jl")

####################################################
# Helper Functions
####################################################

⊗(A,B) = kron(A,B) #just makes things look nice

#Kronecker delta function
#δ(x::T, y::T) where {T<:Number} = isequal(x,y) ? one(T) : zero(T) #this typing is excessive, but good practice 
δ(x, y) = isequal(x,y) ? 1 : 0


#The following functions help generate adjacency matrices for the underlying network. Circulant arrays are used for periodic boundaries and off-diagonal arrays are used for non-periodic graphs
#See https://stackoverflow.com/a/45958661


#Creates a circulant array. Rows are shifted over one down the array
# ⋅  1  ⋅  ⋅  1
# 1  ⋅  1  ⋅  ⋅
# ⋅  1  ⋅  1  ⋅
# ⋅  ⋅  1  ⋅  1
# 1  ⋅  ⋅  1  ⋅
function circulant!(A,n)

    vec = spzeros(Int, n)
    vec[[2,n]] .= 1 #put values into the 2nd and last entry

    for i in 1:n
        A[i,:] = vec
        circshift!(vec, vec, 1) #shift the row over by 1
    end

    return A
end

#Fills the diagonals above and below main diagonal
# ⋅  1  ⋅  ⋅  ⋅
# 1  ⋅  1  ⋅  ⋅
# ⋅  1  ⋅  1  ⋅
# ⋅  ⋅  1  ⋅  1
# ⋅  ⋅  ⋅  1  ⋅
function offDiags!(A)
    A[diagind(A,1)] .= 1
    A[diagind(A,-1)] .= 1
    return A
end

#Generates a square (n,n) sparse adjacency matrix to convert into a network
function genAdj(n::Int, isPeriodic::Bool)
    A = spzeros(Int, n,n)
    
    if isPeriodic #Does the graph wrap around?
        circulant!(A,n)
    else
        offDiags!(A)
    end

    return A ⊗ I(n) + I(n) ⊗ A
end

#Generates a rectangular (m,n) sparse adjacency matrix to convert into a network
function genAdj(m::Int, n::Int, isPeriodic::Bool)
    Am = spzeros(Int, m, m)
    An = spzeros(Int, n, n)

    if isPeriodic #Does the graph wrap around?
        circulant!(Am,m)
        circulant!(An,n)
    else
        offDiags!(Am)
        offDiags!(An)
    end

    return An ⊗ I(m) + I(n) ⊗ Am
end

#Generates a 3D (l,m,n) sparse adjacency matrix to convert into a network
function genAdj(n::Int, m::Int, l::Int, isPeriodic::Bool)
    Al = spzeros(Int, l,l)
    Am = spzeros(Int, m,m)
    An = spzeros(Int, n,n)
    
    if isPeriodic #Does the graph wrap around?
        circulant!(Al,l)
        circulant!(Am,m)
        circulant!(An,n)
    else
        offDiags!(Al)
        offDiags!(Am)
        offDiags!(An)
    end

    return I(l) ⊗ (Am ⊗ I(n) + I(m) ⊗ An) + Al ⊗ I(m*n) #This took a while to figure out...
end



####################################################
# Position the cells
####################################################

##########output is sparse? ⤵

#Creates a cluster of cells in the middle of the grid
function initializeCells(g, gridSize, cellCountsMedium, cellVolumesMedium)

    #Remove Medium from the cells
    cellCounts = cellCountsMedium[2:end]
    cellVolumes = cellVolumesMedium[2:end]

    #initialize matrix of cell IDs (σ)
    cellMembership = zeros(Int, gridSize...)
    
    #Find the center node of the entire graph
    centerIdx = CartesianIndex(gridSize.÷2)
    nodeIdx = LinearIndices(gridSize)
    centerNode = nodeIdx[centerIdx]

    #Determine how far nodes are from the center
    nodeDis = gdistances(g, centerNode)

    #How many nodes need to be initialized?
    totalNodes = cellCounts ⋅ cellVolumes #dot product to sum and multiply

    #Get a sorted permutation of the distance
    sortedDis = sortperm(nodeDis)

    #Assign 1:totalNodes to be filled with cells and the rest medium
    networkIdx = sortedDis[1:totalNodes] 

    #Partition the identified nodes by the number of cells needed
    if sum(cellCounts) == 1 #There is only one cell (no need to partition)
        cellMembership[networkIdx] .= 1
    else
        cellMembership[networkIdx] = Metis.partition(g[networkIdx], sum(cellCounts))
    end

    return cellMembership
end


####################################################
# Model Structure
####################################################

#Each cell in the model is given some basic information like location and size
mutable struct cell
    name::Symbol
    id::Int
    volume::Int
end



#Main container for the model
mutable struct cellpottmodel{N,T}
    gridSize::NTuple{N, Int}
    cells::Vector{cell}
    graph::network{T}
    typeMap::Dict{Symbol, Int}
    energy::Int
    
    function cellpottmodel(;
        gridSize=(50,50),
        cellCounts=[10,5],
        cellVolumes=[20,35],
        cellTypes=[:cell1, :cell2])
        
        #Define the default penalties
        penalties = [AdhesionPenalty(cellTypes),
                     VolumePenalty(cellVolumes, ones(Int,length(cellVolumes)))]
        
        #Add Medium to the cell information
        pushfirst!(cellCounts, 0)
        pushfirst!(cellVolumes, 0)
        pushfirst!(cellTypes, :Medium)


        #How many cells total?
        totalCells = sum(cellCounts)
        
        #Create the data structure for the network
        adjMatrix = genAdj(gridSize..., true)
        g = network(adjMatrix)

        #Initialize cell IDs for each node
        cellIDArray = initializeCells(g, gridSize, cellCounts, cellVolumes)
        
        
        #Create the data structure for the cells (exclude Medium)
        cells = cell.(
            inverse_rle(cellTypes, cellCounts), #inverse_rle(["a","b"], [2,3]) = ["a","a","b","b","b"]
            1:totalCells,
            inverse_rle(cellVolumes, cellCounts),
            [Int[] for _ in 1:totalCells])


        #Update the cells locations, node IDs, and node types
        for (i, cellID) in enumerate(cellIDArray)

            if cellID ≠ 0
                push!(cells[cellID].locations, i)
            end

            push!(g.nodeIDs, cellID)
            push!(g.nodeTypes, cellID == 0 ? :Medium : cells[cellID].name)
        end

        #Associate a number with every cell type (e.g. :Medium→0, CellA→1, CellB→2,...)
        typeMap = Dict(cellTypes .=> 0:length(cellTypes)-1)

        energy = 0

        #Create a new instance of the model
        N = length(gridSize)
        T = eltype(gridSize)
        cpm = new{N,T}(gridSize, cells, g, typeMap, energy)

        #Update the energy penalty
        for penalty in penalties
            totalPenalty!(cpm, penalty)
        end

        return cpm
    end
end


####################################################
# Penalties
####################################################

mutable struct AdhesionPenalty
    J::OffsetMatrix{Int, Matrix{Int}} #J[n,m] gives the adhesion penality for cells with types n and m
    

    function AdhesionPenalty(J::Matrix{Int})
        issymmetric(J) ? nothing : error("J needs to be symmetric")
        
        Joff = OffsetArray(J, CartesianIndex(0, 0):CartesianIndex(size(J).-1))
        return new(Joff)
    end
end

#Just given a list of cell types
function AdhesionPenalty(cellTypes::Vector{Symbol})

    #Count number of cell types including Medium
    uniqueCellCount = length(cellTypes)
    
    #Default adhesion
        # Medium ↔ Medium = 0
        # Medium ↔ Cell = 50
        # Cell ↔ Cell = 100
        # Cell A ↔ Cell B = 75

    J = fill(75, uniqueCellCount, uniqueCellCount)
    J[diagind(J)] .= 100
    J[1,2:end] .= 50
    J[2:end,1] .= 50
    J[1,1] = 0

    return AdhesionPenalty(J)
end


mutable struct VolumePenalty
    desiredVolumes::OffsetVector{Int,Vector{Int64}}
    λᵥ::OffsetVector{Int,Vector{Int64}}

    function VolumePenalty(desiredVolumes::Vector{Int}, λᵥ::Vector{Int})
        desiredVolumesOff = OffsetVector([0; desiredVolumes], 0:length(desiredVolumes))
        λᵥOff = OffsetVector([0; λᵥ], 0:length(λᵥ))
        return new(desiredVolumesOff, λᵥOff)
    end
end


####################################################
# Penalty Functions Updaters
####################################################


function totalPenalty!(cpm::cellpottmodel, AP::AdhesionPenalty)

    #Extract the graph from the CellularPotts Model
    g = cpm.graph

    #Loop through all edges in network and calculate adhesion adjacencies
    for edge in edges(g)

        nodeᵢⱼ = [src(edge), dst(edge)]
        
        #Given a node index, get the cell id and type
        (σᵢ, σⱼ) = g.nodeIDs[nodeᵢⱼ]
        (typeᵢ, typeⱼ) = g.nodeTypes[nodeᵢⱼ]

        #Convert the type (string) to an index for J
        (τᵢ, τⱼ) = cpm.typeMap[typeᵢ], cpm.typeMap[typeⱼ]

        #Only add energy if the nodes are different
        cpm.energy += AP.J[τᵢ, τⱼ] * (1-δ(σᵢ, σⱼ))
    end

    #return H
    return cpm
end


function totalPenalty!(cpm::cellpottmodel, VP::VolumePenalty)

    #Loop through cells, see how far they are from a desired volume, add to energy penality
    for currCell in cpm.cells
        #take the cell's type and find its specific volume constrant (λ)
        typeIdx = cpm.typeMap[currCell.name]
        λ = VP.λᵥ[typeIdx]
        cpm.energy += λ * (currCell.volume - VP.desiredVolumes[typeIdx])^2
    end

    return cpm
end


####################################################
# Override Base.show for each struct
####################################################

function Base.show(io::IO, cpm::cellpottmodel) 
    println("Cell Potts Model:")
    #Grid
    dim = length(cpm.gridSize)
    if dim == 1
        println("$(cpm.gridSize[1])×$(cpm.gridSize[1]) Grid")
    elseif dim == 2
        println("$(cpm.gridSize[1])×$(cpm.gridSize[2]) Grid")
    else
        println("$(cpm.gridSize[1])×$(cpm.gridSize[2])×$(cpm.gridSize[3]) Grid")
    end

    #Cells and types
    cellCounts = countmap([c.name for c in cpm.cells])
    print("Cell Counts:")
    for (key, value) in cellCounts #remove medium
        print(" [$(key) → $(value)]")
    end

    if length(cellCounts) > 1
        println(" [Total → $(length(cpm.cells))]")
    else
        print("\n")
    end

    print("Model Penalties:")
    for p in typeof.(cpm.M.penalties)
        print(" $(p)")
    end
    print("\n")
    println("Current Energy: ", cpm.energy)
    # println("Grid Temperature: ", cpm.M.temperature)
    # println("Steps: ", cpm.stepCounter)
end

