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

#Returns a zero indexed array
offset(x) = OffsetVector(x, Origin(0))

#Check if any lengths in a collection are different
#Works even if the you can't broadcast a function to the collection (e.g. Dict, NamedTuple)
function differentSizes(x)
    validLength = length(first(x))

    for obj in x
        if length(obj) ≠ validLength
            return true
        end
    end

    return false
end

####################################################
# Table to hold cell information
####################################################

#Could use something like TypedDelegation.jl to pass functions to data

mutable struct cellTable{T<:NamedTuple}
    data::T
end

####################################################
# Same basic methods for cellTable
####################################################

countcells(df::cellTable) = length(df.cellIDs) - 1
countcelltypes(df::cellTable) = length(unique(df.typeIDs)) - 1

#From Base, used for stuff like view(Array) to get the Array back
parent(df::cellTable) = getfield(df,:data)

getproperty(df::cellTable, name::Symbol) = getproperty(parent(df), name)

merge(df::cellTable, newColumn) = cellTable( merge(parent(df), newColumn) )

keys(df::cellTable) = keys(parent(df))
values(df::cellTable) = values(parent(df))
pairs(df::cellTable) = pairs(parent(df))

#TODO maybe impliment a row object, or think about Tables.jl more
iterate(df::cellTable, iter=1) = iter ≥ countcells(df) ? nothing : ((;((k,v[iter]) for (k,v) in pairs(df))...), iter + 1)
getindex(df::cellTable, i::Int) = (;((k,[v[i]]) for (k,v) in pairs(df))...)


####################################################
# Function to create a new cell state
####################################################

function newCellState(names::Vector{Symbol}, volumes::Vector{T}, counts::Vector{T}) where T<:Integer

    #Does not include Medium
    totalCells = sum(counts)

    #Add Medium
    pushfirst!(names, :Medium)
    pushfirst!(counts, one(T))
    pushfirst!(volumes, zero(T))

    data =  (;
        names = inverse_rle(names, counts),  #inverse_rle(["a","b"], [2,3]) = ["a","a","b","b","b"] 
        cellIDs = collect(0:totalCells), 
        typeIDs = inverse_rle(0:length(names)-1, counts), 
        volumes = zeros(T,totalCells + 1),
        desiredVolumes = inverse_rle(volumes, counts),
        perimeters = zeros(T,totalCells + 1),
        desiredPerimeters = estPerimeter.(inverse_rle(volumes, counts))
    )

    return cellTable(map(offset, data))
     
end

#Add property for one cell type
function addCellProperty(df::cellTable, propertyName::Symbol, defaultValue, cellName::Symbol)

    newColumn = offset([name == cellName ? defaultValue : missing for name in df.names])

    return merge(df, [propertyName => newColumn])
end

#Or for more than one cell type
function addCellProperty(df::cellTable, propertyName::Symbol, defaultValue, cellName::Vector{Symbol})

    newColumn = offset([name ∈ cellName ? defaultValue : missing for name in df.names])

    return merge(df, [propertyName => newColumn])
end

#Or for all cells
function addCellProperty(df::cellTable, propertyName::Symbol, values::Vector{T}) where T

    pushfirst!(values, first(values)) # note sure what to do about medium
    newColumn = offset(values)

    return merge(df, [propertyName => newColumn])
end


function addNewCell(df::cellTable, cell::T) where T<:NamedTuple
    #TODO Need some checks (e.g. all cells have unique ids, all the keys match)
    for property in keys(df)
        append!(parent(df)[property], cell[property])
    end
end

####################################################
# Variables for Markov Step 
####################################################

#These could be static arrays?
mutable struct MHStepInfo{T<:Integer}
    sourceNode::T                 #Index of node choosen
    neighborNodes::Vector{T}      #Indicies for the neighboring nodes
    sourceCellID::T               #ID of sourceNode
    targetCellID::T               #ID of chosen cell target
    stepCounter::T                #Counts the number of MHSteps performed

    function MHStepInfo(T::DataType)
        return new{T}(zero(T), zeros(T,8), zero(T), zero(T), zero(T))
    end
end


####################################################
# Structure for the model
####################################################

mutable struct CellPotts{N, T<:Integer, V}
    space::CellSpace{N,T}
    initialState::cellTable{V}
    currentState::cellTable{V}
    penalties::Vector{Penalty}
    step::MHStepInfo{T}
    visual::Array{Int,N}
    temperature::Float64

    function CellPotts(space::CellSpace{N,T}, initialCellState::cellTable{V}, penalties::Vector{Penalty}) where {N,T,V}

        return new{N,T,V}(
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

countcells(cpm::CellPotts) = countcells(cpm.currentState)
countcelltypes(cpm::CellPotts) = countcelltypes(cpm.currentState)


####################################################
# Penalty Functors
####################################################

function (AP::AdhesionPenalty)(cpm::CellPotts)
    #calculate Δadhesion to see if the target will decrease model energy
    return AP(cpm, cpm.step.targetCellID) - AP(cpm, cpm.step.sourceCellID)
end

#Moved out of main function to avoid code repeat
function (AP::AdhesionPenalty)(cpm::CellPotts, σᵢ::T) where T<:Integer

    #Initialize the penality
    adhesion = 0

    τᵢ = cpm.currentState.typeIDs[σᵢ]

    for neighbor in cpm.step.neighborNodes
        #Given a node index, get the cellID
        σⱼ = cpm.space.nodeIDs[neighbor]

        #Convert the cellID to cellType
        τⱼ = cpm.currentState.typeIDs[σⱼ]

        #Adhesion is increased if adjacent cells are different types
        adhesion += AP.J[τᵢ, τⱼ] * (1-δ(σᵢ, σⱼ))
    end

    return adhesion
end


function (VP::VolumePenalty)(cpm::CellPotts)

    σᵢ = cpm.step.sourceCellID
    σⱼ = cpm.step.targetCellID
    sourceVolume = VP(cpm, σᵢ) + VP(cpm, σⱼ)
   
    #Change the volumes and recalculate penalty
    cpm.currentState.volumes[σᵢ] -= 1
    cpm.currentState.volumes[σⱼ] += 1

    targetVolume = VP(cpm, σᵢ) + VP(cpm, σⱼ)

    #Reset the volumes
    cpm.currentState.volumes[σᵢ] += 1
    cpm.currentState.volumes[σⱼ] -= 1

    return targetVolume - sourceVolume
end

#Moved out of main function to avoid code repeat
function (VP::VolumePenalty)(cpm::CellPotts, σ::T) where T<:Integer

    volume = cpm.currentState.volumes[σ]
    desiredVolume = cpm.currentState.desiredVolumes[σ]
    τⱼ = cpm.currentState.typeIDs[σ]

    return VP.λᵥ[τⱼ] * (volume - desiredVolume)^2
end

#TODO Add perimeter constraint
function (PP::PerimeterPenalty)(cpm::CellPotts)

    σᵢ = cpm.step.sourceCellID
    σⱼ = cpm.step.targetCellID
    sourcePerimeter = PP(cpm, σᵢ) + PP(cpm, σⱼ)
   
    #Change the volumes and recalculate penalty
    cpm.currentState.volumes[σᵢ] -= 1
    cpm.currentState.volumes[σⱼ] += 1

    targetPerimeter = PP(cpm, σᵢ) + PP(cpm, σⱼ)

    #Reset the volumes
    cpm.currentState.volumes[σᵢ] += 1
    cpm.currentState.volumes[σⱼ] -= 1

    return targetPerimeter - sourcePerimeter
end

#Moved out of main function to avoid code repeat
function (VP::VolumePenalty)(cpm::CellPotts, σ::T) where T<:Integer

    perimeter = cpm.currentState.perimeters[σ]
    desiredPerimeter = cpm.currentState.desiredPerimeters[σ]
    τⱼ = cpm.currentState.typeIDs[σ]

    return VP.λᵥ[τⱼ] * (perimeter - desiredPerimeter)^2
end

####################################################
# Override Base.show
####################################################

function show(io::IO, cpm::CellPotts) 
    println("Cell Potts Model:")
    #Grid
    dim = length(cpm.space.gridSize)
    if dim == 2
        println("Grid: $(cpm.space.gridSize[1])×$(cpm.space.gridSize[2])")
    else
        println("Grid: $(cpm.space.gridSize[1])×$(cpm.space.gridSize[2])×$(cpm.space.gridSize[3])")
    end

    #Cells and types
    cellCounts = countmap(cpm.currentState.names)
    print("Cell Counts:")
    for (key, value) in cellCounts #remove medium
        if key ≠ :Medium
            print(" [$(key) → $(value)]")
        end
    end

    if length(cellCounts) > 1
        println(" [Total → $(length(cpm.currentState.names)-1)]")
    else
        print("\n")
    end

    print("Model Penalties:")
    for p in typeof.(cpm.penalties)
        print(" $(p)")
    end
    print("\n")
    println("Temperature: ", cpm.temperature)
    print("Steps: ", cpm.step.stepCounter)
end

function show(io::IO, intState::cellTable) 

    data = map(OffsetArrays.no_offset_view, getfield(intState,:data))

    hl = Highlighter(f = (data, i, j) -> i == 1,
                         crayon = Crayon(background = :dark_gray))

    pretty_table(
        data,
        header_crayon = crayon"yellow bold",
        highlighters = hl,
        display_size = (20,0),
        vcrop_mode = :middle)
end