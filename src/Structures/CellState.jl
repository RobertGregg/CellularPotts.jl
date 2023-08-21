####################################################
# Table to hold cell information
####################################################

#This is basically a fancy table that has a zero based index
#We want the first (0th) row to correspond with :Medium, this lets cell 1 have index 1 etc.

#Subtyping Tables to 
    #(1) integrate with other like CSV.jl
    #(2) enable fast row/column iteration

"""
    CellState
An concrete type that stores a table where each row is a cell and each column is a cell property.
"""
mutable struct CellState{T<:NamedTuple} <: Tables.AbstractColumns
    data::T
end


"""
    CellRow
An concrete type that stores one row of a CellState.
"""
struct CellRow{T} <: AbstractRow
    row::Int
    source::CellState{T}
end

####################################################
# Tables.jl specific methods
####################################################

# declare that CellState is a table
istable(::Type{<:CellState}) = true


#From Base, used for stuff like view(Array) to get the Array back
parent(df::CellState) = getfield(df,:data)


# schema is column names and types
schema(df::CellState{T}) where T = Schema(keys(parent(df)), eltype.(values(parent(df))) )

# column interface
columnaccess(::Type{<:CellState}) = true
columns(df::CellState) = df
# required AbstractColumns object methods
getcolumn(df::CellState, ::Type{T}, col::Int, nm::Symbol) where {T} = parent(df)[col]
getcolumn(df::CellState, nm::Symbol) = getfield(parent(df), nm)
getcolumn(df::CellState, i::Int) = parent(df)[i]
columnnames(df::CellState) = keys(parent(df))


# required `AbstractRow` interface methods (same as for `AbstractColumns` object before)
# but this time, on our custom row type
getcolumn(r::CellRow, ::Type, col::Int, nm::Symbol) = parent(getfield(r, :source))[col][getfield(r, :row)]
getcolumn(r::CellRow, i::Int) = parent(getfield(r, :source))[i][getfield(r, :row)]
getcolumn(r::CellRow, nm::Symbol) = getfield(parent(getfield(r, :source)), nm)[getfield(r, :row)]
columnnames(df::CellRow) = keys(parent(getfield(df, :source)))


# declare that any CellState defines its own `rows` method
rowaccess(::Type{<:CellState}) = true
# just return itself, which means CellState must iterate `AbstractRow`-compatible objects
rows(df::CellState) = df

####################################################
# Overloading Base Methods
####################################################

eltype(df::CellState{T}) where {T} = CellRow{T}
length(df::CellState) = length(first(parent(df)))
size(df::CellState) = (length(df),length(parent(df)))
iterate(df::CellState, st=0) = st > length(df)-1 ? nothing : (CellRow(st, df), st + 1)
getindex(df::CellState,i::Int) = CellRow(i, df)
merge(df::CellState, newColumn) = CellState( merge(parent(df), newColumn) )


####################################################
# Function to create a new cell state
####################################################

"""
    CellState(names::Vector{Symbol}, volumes::Vector{T}, counts::Vector{T}) where T<:Integer
Create a new `CellState` where each row corresponds to a cell.

By default, this function generates a table with the following columns:
 - names`::Vector{Symbol}` -- List of names given to cells (e.g. `:TCell`)
 - cellIDs`::Vector{<:Integer}` -- A unqiue number given to a cell
 - typeIDs`::Vector{<:Integer}` -- A number corresponding to the cell's name
 - volumes`::Vector{<:Integer}` -- Number of grid squares occupied 
 - desiredVolumes`::Vector{<:Integer}`-- Desired number of grid square
 - perimeters`::Vector{<:Integer}`-- Cell border penality
 - desiredPerimeters`::Vector{<:Integer}`-- Desired cell border penality

The first row in the table is reserved for `:Medium` which is the name given to grid locations not belonging to any cell and is given an index of 0 (The first cell is given an index of 1).
    
Of note, `desiredPerimeters` are calculated as the minimal perimeter given the cell's volume. 
"""
function CellState(names::AbstractVector{S}, volumes::AbstractVector{T}, counts::AbstractVector{T}) where {S<:Union{Symbol, String}, T<:Integer} 

    #Does not include Medium
    totalCells = sum(counts)

    #Add Medium
    pushfirst!(names, first(names) isa Symbol ? :Medium : "Medium")
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

    return CellState(map(offset, data))
     
end

#Alternative method when only creating one cell
CellState(names::Union{Symbol, String}, volumes::T, counts::T) where T<:Integer = CellState([names], [volumes], [counts])

####################################################
# Add/remove cells and properties
####################################################

"""
    addcellproperty(df::CellState, propertyName, propertyValue)
    addcellproperty(df::CellState, propertyName, propertyValue, validCells)
    addcellproperty(df::CellState, propertyName, cellPropertyPairs::Dict{Symbol, T})

Given a `CellState`, add a new column called `propertyName` with values from `propertyValue`.

There are several ways to add a new property to a `CellState`, with the simplest method requiring only a `propertyName` and a single `propertyValue`. A vector of `propertyValue`s can also be provided if a unique value for each cell is desired.

Sometimes a property may only apply to certain cell types. A single cell type or vector of cell types can be passed to `validCells`. All non-valid cells will be given a `missing` value for that property.

Finally, a dictionary of cell=>value pairs can be passed for a given property. Cells not in this dictionary will be given a property value of `missing`.
"""
function addcellproperty(df::CellState, propertyName::Symbol, propertyValue::Vector{T}) where T
    
    #Check if :Medium is accounted for
    if length(propertyValue) ≠ length(df)
        pushfirst!(propertyValue, first(propertyValue))
    end
    propertyValue = offset(propertyValue)

    return merge(df, [propertyName => propertyValue])
end

#Just given and name and value, direct to default method above
addcellproperty(df::CellState, propertyName::Symbol, propertyValue) = addcellproperty(df, propertyName, fill(propertyValue, length(df)))


function addcellproperty(df::CellState, propertyName::Symbol, cellPropertyDict::Dict{Symbol, T}) where T

    newColumn = Union{Missing, T}[]

    for cellName in df.names
        if !haskey(cellPropertyDict, cellName)
            push!(newColumn, missing)
        else
            push!(newColumn, cellPropertyDict[cellName])
        end
    end

    if all(!ismissing, newColumn)
        newColumn = convert(Vector{T}, newColumn)
    end
    
    # direct to default method above
    return addcellproperty(df, propertyName, newColumn)
end

addcellproperty(df::CellState, propertyName::Symbol, propertyValue, validCells) = addcellproperty(df, propertyName, Dict(validCells .=> propertyValue))

function addcellproperty(df::CellState, propertyName::Symbol, propertyValue::Vector{T}, validCells::Symbol) where T
    if length(propertyValue) > 1
        throw(TypeError(Symbol("propertyValue because multiple values cannot apply to one cell"), T, propertyValue))
    end

    return addcellproperty(df, propertyName, Dict(validCells .=> propertyValue))
end



"""
    addnewcell(df::CellState, cell<:NamedTuple)

Given a `CellState`, add a new row corresponding to a new cell in the model. Property names in the for the cell need to match column names in the CellState
"""
function addnewcell(df::CellState, cell::CellRow)
    #TODO Need some checks (e.g. all the keys match)
    for property in keys(df)
        push!(df[property], cell[property])
    end
end

"""
    removecell(df::CellState, cellID)

Given a `CellState`, remove the cell with provided `cellID`.
"""
function removecell(df::CellState, cellID::T) where T<:Integer
    for property in keys(df)
        deleteat!(df[property], cellID)
    end
end


####################################################
# Misc functions
####################################################

countcells(df::CellState) = length(df.cellIDs) - 1
countcelltypes(df::CellState) = length(unique(df.typeIDs)) - 1

####################################################
# Override Base.show
####################################################

#TODO: Add a compact mode

function show(io::IO, intState::CellState) 

    data = map(OffsetArrays.no_offset_view, parent(intState))

    hl = Highlighter(f = (data, i, j) -> i == 1,
                         crayon = Crayon(background = :dark_gray))

    pretty_table(io,
        data,
        header_crayon = crayon"yellow bold",
        highlighters = hl,
        #display_size = (20,0),
        vcrop_mode = :middle)
end
