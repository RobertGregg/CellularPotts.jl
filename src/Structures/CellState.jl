####################################################
# Table to hold cell information
####################################################

#This is basically a fancy table that has a zero based index
#We want the first (0th) row to correspond with :Medium, this lets cell 1 have index 1 etc.

"""
    CellState
An concrete type that stores a table where each row is a cell and each column is a cell property.
"""
struct CellState{T<:NamedTuple}
    data::T
end

####################################################
# Create a row structure for iteration
####################################################

"""
    Cell
An concrete type that points to a row in a referenced CellState table.
"""
struct Cell{T}
    row::Int
    cs::CellState{T}
end

getrow(cell::Cell) = getfield(cell,:row)

Base.parent(cell::Cell) = parent(getfield(cell,:cs))
Base.getproperty(cell::Cell, property::Symbol) = getproperty(parent(cell), property)[getrow(cell)]
Base.propertynames(cell::Cell) = propertynames(parent(cell)) # for tab autocomplete

function Base.show(io::IO, cell::Cell)
    print(io, map(x->x[getrow(cell)], parent(cell)))
end


####################################################
# Overloading Base Methods for CellState
####################################################

Base.parent(cs::CellState) = getfield(cs,:data)
Base.getproperty(cs::CellState, property::Symbol) = getproperty(parent(cs), property)
Base.length(cs::CellState) = length(parent(cs)) #length will be number of columns
Base.size(cs::CellState) = (length(first(parent(cs))), length(cs))
Base.size(cs::CellState, d::Integer) = size(cs)[d]
Base.propertynames(cs::CellState) = propertynames(parent(cs)) # for tab autocomplete

Base.getindex(cs::CellState, i::Int) = Cell(i, cs)
Base.iterate(cs::CellState, row=1) = row > first(size(cs))-1 ? nothing : (Cell(row, cs), row + 1)

#For adding/removing cells
function Base.push!(cs::CellState, cell::Cell)

    for property in propertynames(cs)
        push!(getproperty(cs, property), getproperty(cell, property))
    end

    cs.cellIDs[end] = maximum(cs.cellIDs) + 1

    return nothing
end

function Base.deleteat!(cs::CellState, i)

     map(x->deleteat!(x,i), parent(cs))

     return nothing
end

function Base.show(io::IO, cs::CellState) 

    #pretty_table doesn't like offset vectors in a NamedTuple
    #convert data into an offset array
    data = map(OffsetArrays.no_offset_view, parent(cs))
    mat = offset(reshape(reduce(vcat,values(data)), size(cs)))

    #gray out medium row
    hl = Highlighter(f = (data, i, j) -> i == 0,
        crayon = Crayon(background = :dark_gray))


    pretty_table(io, mat; 
        highlighters = hl,
        show_row_number=true,
        header=(collect(keys(data)), eltype.(values(data))))     
end


####################################################
# Function to create a new cell state
####################################################

"""
    CellState(
    names::AbstractVector{S},
    volumes::AbstractVector{T},
    counts::AbstractVector{T};
    options...) where {S<:Symbol, T<:Integer}


    CellState(;names::AbstractVector{Symbol}, volumes::AbstractVector{T}, counts::AbstractVector{T}, options...) where T<:Integer
Create a new `CellState` where each row corresponds to a cell.

By default, this function generates a table with the following columns:
 - names`::Vector{Symbol}` -- List of names given to cells (e.g. `:TCell`)
 - cellIDs`::Vector{<:Integer}` -- A unqiue number given to a cell
 - typeIDs`::Vector{<:Integer}` -- A number corresponding to the cell's name
 - volumes`::Vector{<:Integer}` -- Number of grid squares occupied 
 - desiredVolumes`::Vector{<:Integer}`-- Desired number of grid squares
 - perimeters`::Vector{<:Integer}`-- Cell border penality
 - desiredPerimeters`::Vector{<:Integer}`-- Desired cell border penality

The first row in the table is reserved for `:Medium` which is the name given to grid locations not belonging to any cell and is given an index of 0 (The first cell is given an index of 1).

Additional cell properties can be supplied as keyword arguements. The length of the keyword arguement needs to match the number of cell types or the total cell count.      
"""
function CellState(
    names::AbstractVector{S},
    volumes::AbstractVector{T},
    counts::AbstractVector{T};
    options...) where {S<:Symbol, T<:Integer} 

    totalCells = sum(counts)
    totalTypes = length(names)

    #Check if all options are correct length
    for (key, value) in options
        if length(value) ≠ totalTypes && length(value) ≠ totalCells
            error("$key is not the correct length, should match number of cells ($totalCells) or cell types ($totalTypes). Got length $(length(value))")
        end
    end

    #Add Medium
    pushfirst!(names, :Medium)
    pushfirst!(counts, one(T))
    pushfirst!(volumes, zero(T))

    #add medium and decode options if needed
    map(values(options)) do option

        pushfirst!(option, first(option))

        if length(option) == length(counts)
            extendedOption = decode(option,counts)
            empty!(option)
            append!(option, extendedOption)
        end
    end

    data =  (;
        names = decode(names, counts),   
        cellIDs = collect(0:totalCells),
        typeIDs = decode(0:totalTypes, counts),
        volumes = zeros(T, totalCells + 1),
        desiredVolumes = decode(volumes, counts),
        perimeters = zeros(T, totalCells + 1),
        desiredPerimeters = estPerimeter.(decode(volumes, counts)),
        options...
    )

    state = CellState(map(offset,data))

    return state
end

#Alternative method when only creating one cell
#options can be one value like :color="red" or a vector with length equal to counts
function CellState(names::Symbol, volumes::T, counts::T; options...) where T<:Integer

    #options is type Base.Pairs and values(options) is a NamedTuple
    optionsVectorized = (; zip( keys(options), vcat.(values(values(options))) )...)

    return CellState([names], [volumes], [counts]; optionsVectorized...)
end

#This is what the user should be calling
CellState(; names, volumes, counts, options...) = CellState(names, volumes, counts; options...)