####################################################
# Helper Functions
####################################################
#Kronecker delta function
δ(x, y) = isequal(x,y) ? 1 : 0

#Given a desired cell volume, calculate minimum perimeter on a square lattice
#There are 3 pages of notes behind this equation
#It's related to the minimum perimeter for a polyomino which is 2ceil(2√V)
#TODO Currently only work for 2D grids
estPerimeter(V::Int) = iszero(V) ? 0 : 4ceil(Int,2sqrt(V)-3) + 2ceil(Int,2sqrt(V+1)-4) + 14



####################################################
# Function to create a new cell state
####################################################

function newCellState(names::Vector{Symbol}, volumes::Vector{T}, counts::Vector{T}) where T<:Integer

    totalCells = sum(counts)

    return DataFrame(
        names = inverse_rle(names, counts),
        ids = 1:totalCells,
        volumes = zeros(T,totalCells),
        desiredVolumes = inverse_rle(volumes, counts),
        perimeters = zeros(T,totalCells),
        desiredPerimeters = estPerimeter.(inverse_rle(volumes, counts))
    )
end

#Add property for one cell type
function addCellProperty!(df::DataFrame, propertyname::Symbol, defaultValue, cellName::Symbol)

    df[!,propertyname] = [name == cellName ? defaultValue : missing for name in df.names]
end

#Or for more than one cell type
function addCellProperty!(df::DataFrame, propertyname::Symbol, defaultValue, cellName::Vector{Symbol})

    df[!,propertyname] = [name ∈ cellName ? defaultValue : missing for name in df.names]
end

####################################################
# Variables for Markov Step 
####################################################

#These could be static arrays?
mutable struct MHStepInfo{T<:Integer}
    sourceNode::T               #Index of node choosen
    sourceNeighbors::Vector{T}  #Indicies for the neighboring nodes
    sourceCell::T               #ID of sourceNode
    targetCell::T               #ID of chosen cell target
    sourceTargetCell::Vector{T} #Combine the source and target together
    stepCounter::T              #Counts the number of MHSteps performed

    function MHStepInfo(T::DataType)
        return new{T}(zero(T), zeros(T,8), zero(T), zero(T), zeros(T,2), zero(T))
    end
end


####################################################
# Structure for the model
####################################################

mutable struct CellPotts{N, T<:Integer}
    space::CellSpace{N,T}
    initialState::DataFrame
    currentState::DataFrame
    penalties::Vector{Penalty}
    step::MHStepInfo{T}
    visual::Array{Int,N}
    temperature::Float64

    function CellPotts(space::CellSpace{N,T}, initialCellState::DataFrame, penalties::Vector{Penalty}) where {N,T}

        return new{N,T}(
            space,
            initialCellState,
            initialCellState,
            penalties,
            MHStepInfo(T),
            zeros(T,space.gridSize),
            20.0)
    end
end

####################################################
# Helper functions for CellPotts
####################################################

countCells(cpm::CellPotts) = nrow(cpm.currentState)
countCellTypes(cpm::CellPotts) = length(unique(cpm.currentState.names))

####################################################
# Penalty Functors
####################################################

#TODO Do we really need OffSetArrays?
function (AP::AdhesionPenalty)(cpm::CellPotts)
    step = cpm.step

    sourceAdhesion = 0
    σᵢ = step.sourceCell
    typeᵢ = cpm.currentState.names[σᵢ]

    for sourceNodeNeighbor in step.sourceNeighbors
        #Given a node index, get the cell id and type
        σⱼ = cpm.space.nodeIDs[sourceNodeNeighbor]
        typeⱼ = cpm.space.nodeTypes[sourceNodeNeighbor]

        #Convert the type (string) to an index for J
        (τᵢ, τⱼ) = ( cpm.cells.typeMap[typeᵢ], cpm.cells.typeMap[typeⱼ] )

        sourceAdhesion += AP.J[τᵢ, τⱼ] * (1-δ(σᵢ, σⱼ))
    end
end

####################################################
# Override Base.show cpm
####################################################

function Base.show(io::IO, cpm::CellPotts) 
    println("Cell Potts Model:")
    #Grid
    dim = length(cpm.space.gridSize)
    if dim == 1
        println("Grid: $(cpm.space.gridSize[1])×$(cpm.space.gridSize[1])")
    elseif dim == 2
        println("Grid: $(cpm.space.gridSize[1])×$(cpm.space.gridSize[2])")
    else
        println("Grid: $(cpm.space.gridSize[1])×$(cpm.space.gridSize[2])×$(cpm.space.gridSize[3])")
    end

    #Cells and types
    cellCounts = countmap(cpm.currentState.names)
    print("Cell Counts:")
    for (key, value) in cellCounts #remove medium
            print(" [$(key) → $(value)]")
    end

    if length(cellCounts) > 1
        println(" [Total → $(length(cpm.currentState.names))]")
    else
        print("\n")
    end

    print("Model Penalties:")
    for p in typeof.(cpm.penalties)
        print(" $(p)")
    end
    print("\n")
    println("Temperature: ", cpm.temperature)
    println("Steps: ", cpm.step.stepCounter)
end