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

####################################################
# Table to hold cell information
####################################################

#Could use something like TypedDelegation.jl to pass functions to data
"""
    CellTable
An concrete type that stores a table where each row is a cell and each column is a cell property.

This table can be genereated using the `newCellState()` function.
"""
mutable struct CellTable{T<:NamedTuple}
    data::T
end

####################################################
# Same basic methods for CellTable
####################################################

countcells(df::CellTable) = length(df.cellIDs) - 1
countcelltypes(df::CellTable) = length(unique(df.typeIDs)) - 1

#From Base, used for stuff like view(Array) to get the Array back
parent(df::CellTable) = getfield(df,:data)

getproperty(df::CellTable, name::Symbol) = getproperty(parent(df), name)

merge(df::CellTable, newColumn) = CellTable( merge(parent(df), newColumn) )

keys(df::CellTable) = keys(parent(df))
values(df::CellTable) = values(parent(df))
pairs(df::CellTable) = pairs(parent(df))

#TODO use Tables.jl
iterate(df::CellTable, iter=1) = iter ≥ countcells(df) ? nothing : ((;((k,v[iter]) for (k,v) in pairs(df))...), iter + 1)
getindex(df::CellTable, i::Int) = (;((k,[v[i]]) for (k,v) in pairs(df))...)


####################################################
# Function to create a new cell state
####################################################

#TODO rename as "createNewTable"
"""
    newCellState(names, volumes, counts)
Create a new `cellTable` where each row corresponds to a cell.

By default, this function generates a table with the following columns:
 - names`::Vector{Symbol}`: List of names given to cells (e.g. `:TCell`)
 - cellIDs`::Vector{<:Integer}`: A unqiue number given to a cell
 - typeIDs`::Vector{<:Integer}`: A number corresponding to the cell's name
 - volumes`::Vector{<:Integer}`: Number of grid squares occupied 
 - desiredVolumes`::Vector{<:Integer}`: Desired number of grid square
 - perimeters`::Vector{<:Integer}`: Cell border penality
 - desiredPerimeters`::Vector{<:Integer}`: Desired cell border penality

The first row in the table is reserved for `:Medium` which is the name given to grid locations not belonging to any cell and is given an index of 0 (The first cell is given an index of 1).
    
Of note, `desiredPerimeters` are calculated as the minimal perimeter given the cell's volume. 
"""
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

    return CellTable(map(offset, data))
     
end

#TODO make the different options for addCellProperty cleaner
#Add property for one cell type
"""
    addCellProperty(df::CellTable, propertyName, defaultValue, cellName)
    addCellProperty(df::CellTable, propertyName, values)

Given a `cellTable`, add a new column called `propertyName` with a given default value.

A `cellTable` is generated using the `newCellState()` function.

`cellName` can be the name of one cell type or a vector of cell types. Cells not included in `cellName` will be given a value of `missing`.

If cellNames are not specified, a vector of `values` can be supplied for every cell type.  
"""
function addCellProperty(df::CellTable, propertyName::Symbol, defaultValue, cellName::Symbol)

    newColumn = offset([name == cellName ? defaultValue : missing for name in df.names])

    return merge(df, [propertyName => newColumn])
end

#Or for more than one cell type
function addCellProperty(df::CellTable, propertyName::Symbol, defaultValue, cellName::Vector{Symbol})

    newColumn = offset([name ∈ cellName ? defaultValue : missing for name in df.names])

    return merge(df, [propertyName => newColumn])
end

#Or for all cells
function addCellProperty(df::CellTable, propertyName::Symbol, values::Vector{T}) where T

    pushfirst!(values, first(values)) # note sure what to do about medium
    newColumn = offset(values)

    return merge(df, [propertyName => newColumn])
end


"""
    addNewCell(df::CellTable, cell<:NamedTuple)

Given a `cellTable`, add a new row corresponding to a new cell in the mode. Property names in the for the cell need to match column names in the cellTable
"""
function addNewCell(df::CellTable, cell::T) where T<:NamedTuple
    #TODO Need some checks (e.g. all cells have unique ids, all the keys match)
    for property in keys(df)
        append!(parent(df)[property], cell[property])
    end
end

"""
    removeCell(df::CellTable, cellID)

Given a `cellTable`, remove the cell with provided `cellID`.
"""
function removeCell(df::CellTable, cellID::T) where T<:Integer
    for property in keys(df)
        deleteat!(parent(df)[property], cellID)
    end
end

####################################################
# Penalties
####################################################

"""
    Penalty
An abstract type representing a constraint imposed onto the cellular potts model.

To add a new penalty, a new struct subtyping `Penalty` needs to be defined and the `addPenalty!()` function needs to be extended to include the new penalty.

**Note**: variables associated with a new penalty may need to be offset such that index 0 maps to :Medium, index 1 maps to :Cell1, etc.
"""
abstract type Penalty end

"""
    AdhesionPenalty(J::Matrix{Int})
An concrete type that penalizes neighboring grid locations from different cells.

Requires a symmetric matrix `J` where `J[n,m]` gives the adhesion penality for cells with types n and m. `J` is zero-indexed meaning `J[0,1]` and `J[1,0]` corresponds to the `:Medium` ↔ `:Cell1` adhesion penalty.

**Note**: `J` is automatically transformed to be a zero-indexed offset array.
"""
struct AdhesionPenalty <: Penalty
    J::OffsetMatrix{Int, Matrix{Int}}

    function AdhesionPenalty(J::Matrix{Int})
        issymmetric(J) ? nothing : error("J needs to be symmetric")
        
        return new(offset(J))
    end
end

"""
    VolumePenalty(λᵥ::Vector{Int})
An concrete type that penalizes cells that deviate from their desired volume.

Requires a vector `λᵥ` with n penalties where n is the number of cell types. `λᵥ` is zero-indexed meaning `λᵥ[0]` corresponds to the `:Medium` volume penalty (which is set to zero).

**Note**: `λᵥ` is automatically transformed to be a zero-indexed offset array and does not require the volume penalty for `:Medium`.
"""
struct VolumePenalty <: Penalty
    λᵥ::OffsetVector{Int,Vector{Int}}

    function VolumePenalty(λᵥ::Vector{Int})
        λᵥOff = offset([0; λᵥ])
        return new(λᵥOff)
    end
end

"""
    PerimeterPenalty(λᵥ::Vector{Int})
An concrete type that penalizes cells that deviate from their desired perimeter.

Requires a vector `λₚ` with n penalties where n is the number of cell types. `λₚ` is zero-indexed meaning `λₚ[0]` corresponds to the `:Medium` perimeter penalty (which is set to zero).

**Note**: `λₚ` is automatically transformed to be a zero-indexed offset array and does not require the perimeter penalty for `:Medium`.
"""
mutable struct PerimeterPenalty <: Penalty
    λₚ::OffsetVector{Int,Vector{Int}}
    Δpᵢ::Int
    Δpⱼ::Int

    function PerimeterPenalty(λₚ::Vector{Int}) 
        λₚOff = offset([0; λₚ])
        return new(λₚOff, 0, 0)
    end
end

#TODO Cells get stuck in other cell's nodeMemory
"""
    MigrationPenalty(maxAct, λ, gridSize)
An concrete type that encourages cells to protude and drag themselves forward.

Two integer parameters control how cells protude:
 - `maxAct`: A maximum activity a grid location can have
 - `λ`: A parameter that controls the strength of this penalty

Increasing `maxAct` will cause grid locations to more likely protrude. Increasing `λ` will cause those protusions to reach farther away. 

`MigrationPenalty` also requires a list of cell types to apply the penalty to and the grid size (Space.gridSize).
"""
mutable struct MigrationPenalty <: Penalty
    maxAct::Int
    λ::OffsetVector{Int,Vector{Int}}
    nodeMemory::SparseMatrixCSC{Int,Int}

    function MigrationPenalty(maxAct::T, λ::Vector{T}, gridSize::NTuple{N,T}) where {T<:Integer, N}
        λOff = offset([0; λ])
        return new(maxAct, λOff, spzeros(T,gridSize))
    end
end



####################################################
# Variables for Markov Step 
####################################################

mutable struct MHStepInfo{T<:Integer}
    sourceNode::T      #Index of node choosen
    targetNode::T      #Index of node choosen
    sourceNeighborNodes::Vector{T} #Indicies for the neighboring nodes
    targetNeighborNodes::Vector{T} #Indicies for the neighboring nodes
    sourceCellID::T    #ID of sourceNode
    targetCellID::T    #ID of chosen cell target
    stepCounter::T     #Counts the number of MHSteps performed
end 

MHStepInfo() = MHStepInfo(0,0,[0],[0],0,0,0)


####################################################
# Structure for the model
####################################################

"""
    CellPotts(space, initialCellState, penalties)
A data container that holds information to run the cellular potts simulation.

Requires three inputs:
 - `space`: a region where cells can exist, generated using `CellSpace()`.
 - `initialCellState`: a table where rows are cells and columns are cell properties, generated using `newCellState()`.
 - `penalties`: a vector of penalties to append to the model.
"""
mutable struct CellPotts{N, T<:Integer, V<:NamedTuple, U}
    space::CellSpace{N,T}
    initialState::CellTable{V}
    currentState::CellTable{V}
    penalties::Vector{U}
    step::MHStepInfo{T}
    getArticulation::ArticulationUtility
    temperature::Float64

    function CellPotts(space::CellSpace{N,T}, initialCellState::CellTable{V}, penalties::Vector{P}) where {N,T,V,P}

        #See https://github.com/JuliaLang/julia/pull/44131 for why Unions are used
        U = Union{typeof.(penalties)...}
        return new{N,T,V,U}(
            space,
            initialCellState,
            initialCellState,
            U[p for p in penalties],
            MHStepInfo(),
            ArticulationUtility(nv(space)),
            20.0)
    end
end

####################################################
# Helper functions for CellPotts
####################################################

"""
    countcells(cpm::CellPotts)
    countcells(df::CellTable)

Count the number of cells in the model 
"""
countcells(cpm::CellPotts) = countcells(cpm.currentState)

"""
    countcelltypes(cpm::CellPotts)
    countcelltypes(df::CellTable)

Count the number of cell types in the model 
"""
countcelltypes(cpm::CellPotts) = countcelltypes(cpm.currentState)

####################################################
# Override Base.show
####################################################

function show(io::IO, cpm::CellPotts) 
    println(io,"Cell Potts Model:")
    #Grid
    dim = length(cpm.space.gridSize)
    if dim == 2
        println(io,"Grid: $(cpm.space.gridSize[1])×$(cpm.space.gridSize[2])")
    else
        println(io,"Grid: $(cpm.space.gridSize[1])×$(cpm.space.gridSize[2])×$(cpm.space.gridSize[3])")
    end

    #Cells and types
    cellCounts = countmap(cpm.currentState.names)
    print(io,"Cell Counts:")
    for (key, value) in cellCounts #remove medium
        if key ≠ :Medium
            print(io," [$(key) → $(value)]")
        end
    end

    if length(cellCounts) > 1
        println(io," [Total → $(length(cpm.currentState.names)-1)]")
    else
        print(io,"\n")
    end

    print(io,"Model Penalties:")
    for p in Base.uniontypes(eltype(cpm.penalties))
        p = Symbol(p)
        print(io," $(replace(String(p),"Penalty"=>""))")
    end
    print(io,"\n")
    println(io,"Temperature: ", cpm.temperature)
    print(io,"Steps: ", cpm.step.stepCounter)
end

function show(io::IO, intState::CellTable) 

    data = map(OffsetArrays.no_offset_view, getfield(intState,:data))

    hl = Highlighter(f = (data, i, j) -> i == 1,
                         crayon = Crayon(background = :dark_gray))

    pretty_table(io,
        data,
        header_crayon = crayon"yellow bold",
        highlighters = hl,
        display_size = (20,0),
        vcrop_mode = :middle)
end