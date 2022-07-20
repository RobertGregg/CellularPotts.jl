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
offset(x) = OffsetArray(x, Origin(0))

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
        typeIDs = inverse_rle(0:length(names)-1, counts), #look into enum?
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
# Penalties
####################################################

abstract type Penalty end

#Offset arrays so the zeroth index refers to Medium

struct AdhesionPenalty <: Penalty
    J::OffsetMatrix{Int, Matrix{Int}} #J[n,m] gives the adhesion penality for cells with types n and m

    function AdhesionPenalty(J::Matrix{Int})
        issymmetric(J) ? nothing : error("J needs to be symmetric")
        
        return new(offset(J))
    end
end

struct VolumePenalty <: Penalty
    λᵥ::OffsetVector{Int,Vector{Int}}

    function VolumePenalty(λᵥ::Vector{Int})
        λᵥOff = offset([0; λᵥ])
        return new(λᵥOff)
    end
end

mutable struct PerimeterPenalty <: Penalty
    λₚ::OffsetVector{Int,Vector{Int}}
    Δpᵢ::Int
    Δpⱼ::Int

    function PerimeterPenalty(λᵥ::Vector{Int}) 
        λₚOff = offset([0; λᵥ])
        return new(λₚOff)
    end
end

mutable struct MigrationPenalty <: Penalty
    maxAct::Int
    λ::Int
    cellTypes::Vector{Symbol}
    nodeMemory::SparseVector{Int,Int}

    function MigrationPenalty(maxAct::T, λ::T, cellTypes::Vector{S}, gridSize::NTuple{N,T}) where {T<:Integer, S<:Symbol, N}
        return new(maxAct, λ, cellTypes, spzeros(T,prod(gridSize)))
    end
end



####################################################
# Variables for Markov Step 
####################################################

#These could be static arrays?
mutable struct MHStepInfo{T<:Integer}
    sourceNode::T                 #Index of node choosen
    targetNode::T                 #Index of node choosen
    neighborNodes::Vector{T}      #Indicies for the neighboring nodes
    sourceCellID::T               #ID of sourceNode
    targetCellID::T               #ID of chosen cell target
    stepCounter::T                #Counts the number of MHSteps performed

    function MHStepInfo(T::DataType)
        return new{T}(zero(T), zero(T), zeros(T,8), zero(T), zero(T), zero(T))
    end
end


####################################################
# Structure for the model
####################################################

mutable struct CellPotts{N, T<:Integer, V<:NamedTuple, P<:Penalty}
    space::CellSpace{N,T}
    initialState::cellTable{V}
    currentState::cellTable{V}
    penalties::Dict{Symbol,P}
    step::MHStepInfo{T}
    visual::Array{Int,N}
    temperature::Float64

    function CellPotts(space::CellSpace{N,T}, initialCellState::cellTable{V}, penalties::Dict{Symbol,P}) where {N,T,V,P}

        return new{N,T,V,P}(
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
    for p in keys(cpm.penalties)
        print(" $(replace(String(p),"Penalty"=>""))")
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