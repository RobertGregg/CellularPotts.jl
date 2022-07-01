include("GraphStructure.jl")

####################################################
# Helper Functions
####################################################
#Kronecker delta function
δ(x, y) = isequal(x,y) ? 1 : 0


####################################################
# Generate grid network
####################################################

#The mod1 function always returns a number between 1 and n
#0 gets mapped to n
#n+1 gets mapped to 1
#This is exactly what is needed for periodic boundaries
#Make mod1 work for CartesianIndex (I is the tuple of indices)
Base.mod1(x::CartesianIndex{N}, i::Int...) where N = CartesianIndex(mod1.(x.I,i))

#Get the 8 neighboring nodes around current position
#Only works for 2D rectangular grids for now
function mooreNeighbors(J::CartesianIndex{N}, gridSize::NTuple{N, Int}) where N
    Δx = CartesianIndex(1,0)
    Δy = CartesianIndex(0,1)
    
    return mod1.([J-Δx-Δy,J-Δy,J+Δx-Δy,
                  J-Δx,        J+Δx,
                  J-Δx+Δy,J+Δy,J+Δx+Δy], gridSize...)
end

#Given a network, add in all the needed edges
function connectGraph!(g::network{Int}, gridSize::NTuple{N, Int}) where N

    grid = reshape(1:prod(gridSize), gridSize...)

    for J in CartesianIndices(grid)
        for i in mooreNeighbors(J,gridSize)
            add_edge!(g, grid[J], grid[i])
        end
    end

    return g
end

####################################################
# Structure to hold info about initial cells
####################################################

#Provide some defaults for an example
Base.@kwdef struct InitialCellState
    name::Vector{Symbol} = [:Epidermal]
    id::AbstractVector{Int} = [1]
    count::Vector{Int} = [100]
    volume::Vector{Int} = [600]
end

#Method where user doesn't need to worry about id
InitialCellState(name::Vector{Symbol}, count::Vector{Int}, volume::Vector{Int}) = InitialCellState(names, 1:length(name), count, volume)

####################################################
# Default Penalties
####################################################
abstract type Penalty end

#Offset arrays so the zeroth index referes to Medium

mutable struct AdhesionPenalty <: Penalty
    J::OffsetMatrix{Int, Matrix{Int}} #J[n,m] gives the adhesion penality for cells with types n and m

    function AdhesionPenalty(J::Matrix{Int})
        issymmetric(J) ? nothing : error("J needs to be symmetric")
        
        Joff = OffsetArray(J, CartesianIndex(0, 0):CartesianIndex(size(J).-1))
        return new(Joff)
    end
end

#
mutable struct VolumePenalty <: Penalty
    λᵥ::OffsetVector{Int,Vector{Int}}

    function VolumePenalty(λᵥ::Vector{Int})
        λᵥOff = OffsetVector([0; λᵥ], 0:length(λᵥ))
        return new(λᵥOff)
    end
end

####################################################
# Structure for user to provide information
####################################################

Base.@kwdef struct Parameters{N} #N is graph dimension (makes types stable)
    gridSize::NTuple{N, Int} = (250,250) #What are the dimensions of the space?
    isPeriodic::Bool = true #Do you want the space to wrap around?
    state::InitialCellState = InitialCellState() #Given some info about the cells
    penalties::Vector{Penalty} = [AdhesionPenalty([0 20; 20 100]), VolumePenalty([50])] #What penalties are added to the model?
    temperature::Float64 = 20.0 #Temperature to weight site flips (larger temp == more acceptances)
end

####################################################
# Summarize all the cells in the model
####################################################

mutable struct CellSummary
    #Vectors are offset to include medium (medium gets index 0, cell 1 gets index 1, etc.)
    names::OffsetVector{Symbol,Vector{Symbol}} #vector of cell type names
    ids::OffsetVector{Int,Vector{Int}} #vector of cell IDs
    volumes::OffsetVector{Int,Vector{Int}} #vector of cell volumes 
    desiredVolumes::OffsetVector{Int,Vector{Int}} #vector of desired cell volumes
    typeMap::Dict{Symbol, Int} #mapping cell types (e.g. "Medium") to an index (e.g. 0)

    function CellSummary(parameters::Parameters)
        #Unpack the initial state
        state = parameters.state

        #How many cells are there?
        totalCells = sum(state.count)

        ids = OffsetVector(collect(0:totalCells), 0:totalCells)

        volumes = OffsetVector(zeros(Int,totalCells+1), 0:totalCells)

        names = OffsetVector([:Medium; inverse_rle(state.name, state.count)], 0:totalCells)

        desiredVolumes = OffsetVector([0; inverse_rle(state.volume, state.count)], 0:totalCells)

        typeMap = Dict( state.name .=> 1:length(state.name) )
        typeMap[:Medium] = 0

        return new(names, ids, volumes, desiredVolumes, typeMap)
    end
end

####################################################
# Variables for Markov Step 
####################################################

#These could be static arrays?
mutable struct MHStepInfo
    sourceNode::Int #Index of node choosen
    sourceNodeNeighbors::Vector{Int} # Indicies for the neighboring nodes
    possibleCellTargets::Vector{Int} # Unique cell IDs of the neighboring nodes
    sourceCell::Int #ID of sourceNode
    targetCell::Int #ID of chosen cell target
    sourceTargetCell::Vector{Int} #Combine the source and target together
end

####################################################
# Structure for the model
####################################################

mutable struct CellPotts{N}
    parameters::Parameters{N}
    graph::network{Int}
    cells::CellSummary
    stepInfo::MHStepInfo
    visual::Array{Int,N} #An array of cell memberships for plotting
    stepCounter::Int #counts the number of MHSteps performed
    energy::Int

    function CellPotts(parameters::Parameters{N}) where N

        graph = network(prod(parameters.gridSize))
        connectGraph!(graph, parameters.gridSize)

        cells = CellSummary(parameters)

        stepInfo = MHStepInfo(1,[1],[1],1,1,[1,1]) #arbitrary default parameters

        visual = zeros(Int, parameters.gridSize)

        return new{N}(parameters, graph, cells, stepInfo, visual, 0, 0)
    end
end

####################################################
# Create cells for network
####################################################

#Creates a cluster of cells in the middle of the grid
function initializeCells!(cpm::CellPotts)

    #Unpack the initial state and parameters
    parameters = cpm.parameters
    state = parameters.state

    #initialize matrix of cell IDs (σ)
    cellMembership = zeros(Int, parameters.gridSize...)
    
    #Find the center node of the entire graph
    centerIdx = CartesianIndex(parameters.gridSize.÷2)
    nodeIdx = LinearIndices(parameters.gridSize)
    centerNode = nodeIdx[centerIdx]

    #Determine how far nodes are from the center
    nodeDis = gdistances(cpm.graph, centerNode)

    #How many nodes need to be initialized?
    totalNodes = state.count ⋅ state.volume #dot product to sum and multiply

    #Get a sorted permutation of the distance
    sortedDis = sortperm(nodeDis)

    #Assign 1:totalNodes to be filled with cells and the rest medium
    networkIdx = sortedDis[1:totalNodes] 

    #Partition the identified nodes by the number of cells needed
    if sum(state.count) == 1 #There is only one cell (no need to partition)
        cellMembership[networkIdx] .= 1
    else
        cellMembership[networkIdx] = Metis.partition(cpm.graph[networkIdx], sum(state.count))
    end

    #Update the network with the new cell locations
    for (i, cellID) in enumerate(cellMembership)
        if cellID ≠ 0
            cpm.graph.nodeIDs[i] = cellID
            cpm.graph.nodeTypes[i] = cpm.cells.names[cellID]

            #Also update the cell summary volumes while you're at it
            cpm.cells.volumes[cellID] += 1
        end
    end


    #Update the intial penalty energy
    for penalty in cpm.parameters.penalties
        penalty(cpm)
    end

    #Fill the the array for the visual
    cpm.visual = cellMembership

    return cpm
end


####################################################
# Penalty Functors
####################################################

#=
Need to define two methods for each penality: 
    - calculate penality on the entire grid
    - calculate change in penality after markov step
=#

#----Calculate on entire Model (for initialization)----


#This code gets repeated a lot so I made it a separate function
function (AP::AdhesionPenalty)(cpm::CellPotts)

    #unpack the cell type → id dictionary
    typeMap = cpm.cells.typeMap

    #Loop through all edges in network and calculate adhesion adjacencies
    for edge in edges(cpm.graph)
        
        #Given a node index, get the cell id and type
        (σᵢ, σⱼ) = cpm.graph.nodeIDs[ [edge.src, edge.dst] ]
        (typeᵢ, typeⱼ) = cpm.graph.nodeTypes[ [edge.src, edge.dst] ]

        #Convert the type (string) to an index for J
        (τᵢ, τⱼ) = ( typeMap[typeᵢ], typeMap[typeⱼ] )

        cpm.energy += AP.J[τᵢ, τⱼ] * (1-δ(σᵢ, σⱼ))
    end

    return cpm
end


function (VP::VolumePenalty)(cpm::CellPotts)
    #unpack cells
    cell = cpm.cells

    #Loop through cells, see how far they are from a desired volume, and put into appropriate slot
    for (name, volume, desiredVolume) in zip(cell.names, cell.volumes, cell.desiredVolumes)
        typeID = cell.typeMap[name]
        cpm.energy += VP.λᵥ[typeID] * (volume - desiredVolume)^2
    end

    #Calculate the penalty
    return cpm
end


#----Calculate after markov step----

#TODO: think of a clever way to avoid code repeat
function (AP::AdhesionPenalty)(cpm::CellPotts, stepInfo::MHStepInfo)

    sourceH = 0
    σᵢ = stepInfo.sourceNode
    typeᵢ = cpm.graph.nodeTypes[σᵢ]
    for sourceNodeNeighbor in stepInfo.sourceNodeNeighbors
        #Given a node index, get the cell id and type
        σⱼ = cpm.graph.nodeIDs[sourceNodeNeighbor]
        typeⱼ = cpm.graph.nodeTypes[sourceNodeNeighbor]

        #Convert the type (string) to an index for J
        (τᵢ, τⱼ) = ( cpm.cells.typeMap[typeᵢ], cpm.cells.typeMap[typeⱼ] )

        sourceH += AP.J[τᵢ, τⱼ] * (1-δ(σᵢ, σⱼ))
    end


    targetH = 0
    σᵢ = stepInfo.targetCell
    typeᵢ = cpm.cells.names[σᵢ] #replace type with new target type
    for sourceNodeNeighbor in stepInfo.sourceNodeNeighbors
        #Given a node index, get the cell id and type
        σⱼ= cpm.graph.nodeIDs[sourceNodeNeighbor]

        typeⱼ = cpm.graph.nodeTypes[sourceNodeNeighbor]

        #Convert the type (string) to an index for J
        (τᵢ, τⱼ) = (cpm.cells.typeMap[typeᵢ], cpm.cells.typeMap[typeⱼ] )

        targetH += AP.J[τᵢ, τⱼ] * (1-δ(σᵢ, σⱼ))
    end

    #Calculate the adhesion difference between target and source
    return targetH - sourceH
end



function (VP::VolumePenalty)(cpm::CellPotts, stepInfo::MHStepInfo)

    sourceTargetCell = stepInfo.sourceTargetCell
    volumes = cpm.cells.volumes[sourceTargetCell]
    
    
    #Current Penalty
    sourceVolume = 0
    
    for (i, σᵢ) in enumerate(sourceTargetCell)
        volume = volumes[i]
        desiredVolume = cpm.cells.desiredVolumes[σᵢ]
        name = cpm.cells.names[σᵢ]
        typeID = cpm.cells.typeMap[name]
        
        sourceVolume += VP.λᵥ[typeID] * (volume - desiredVolume)^2
    end
    
    #New Penalty
    targetVolume = 0
    #Update the volumes
    volumes[1] -= 1 
    volumes[2] += 1 

    for (i, σᵢ) in enumerate(sourceTargetCell)
        volume = volumes[i]
        desiredVolume = cpm.cells.desiredVolumes[σᵢ]
        name = cpm.cells.names[σᵢ]
        typeID = cpm.cells.typeMap[name]

        targetVolume += VP.λᵥ[typeID] * (volume - desiredVolume)^2
    end
  

    return targetVolume - sourceVolume
end


####################################################
# Override Base.show cpm
####################################################

function Base.show(io::IO, cpm::CellPotts) 
    println("Cell Potts Model:")
    #Grid
    dim = length(cpm.parameters.gridSize)
    if dim == 1
        println("$(cpm.parameters.gridSize[1])×$(cpm.parameters.gridSize[1]) Grid")
    elseif dim == 2
        println("$(cpm.parameters.gridSize[1])×$(cpm.parameters.gridSize[2]) Grid")
    else
        println("$(cpm.parameters.gridSize[1])×$(cpm.parameters.gridSize[2])×$(cpm.parameters.gridSize[3]) Grid")
    end

    #Cells and types
    cellCounts = countmap(cpm.cells.names)
    print("Cell Counts:")
    for (key, value) in cellCounts #remove medium
        if key ≠ :Medium
            print(" [$(key) → $(value)]")
        end
    end

    if length(cellCounts) > 1
        println(" [Total → $(length(cpm.cells.names))]")
    else
        print("\n")
    end

    print("Model Penalties:")
    for p in typeof.(cpm.parameters.penalties)
        print(" $(p)")
    end
    print("\n")
    println("Current Energy: ", cpm.energy)
    # println("Grid Temperature: ", cpm.M.temperature)
    # println("Steps: ", cpm.stepCounter)
end