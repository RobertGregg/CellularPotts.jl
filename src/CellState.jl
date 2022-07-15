#This is basically a fancy table that has a zero based index
#We want the first (0th) row to correspond with :Medium, this lets cell 1 have index 1 etc.

#TODO getindex method?

#TODO T leads to abstract typing, need macro?
mutable struct cellTable{T} <: AbstractColumns
    columnNames::Vector{Symbol}
    lookup::Dict{Symbol, Int}
    data::T
end

# declare that cellTable is a table
istable(::Type{<:cellTable}) = true

# schema is column names and types
schema(df::cellTable{T}) where T = Schema(getfield(df, :columnNames), eltype.(getfield(df, :data)))

# column interface
columnaccess(::Type{<:cellTable}) = true
columns(df::cellTable) = df
# required AbstractColumns object methods
getcolumn(df::cellTable, ::Type{T}, col::Int, nm::Symbol) where {T} = getfield(df, :data)[col]
getcolumn(df::cellTable, nm::Symbol) = getfield(df, :data)[getfield(df, :lookup)[nm]]
getcolumn(df::cellTable, i::Int) = getfield(df, :data)[i]
columnnames(df::cellTable) = getfield(df, :columnNames)


# declare that any cellTable defines its own `rows` method
rowaccess(::Type{<:cellTable}) = true
# just return itself, which means cellTable must iterate `AbstractRow`-compatible objects
rows(df::cellTable) = df
# the iteration interface, at a minimum, requires `eltype`, `length`, and `iterate`
# for `cellTable` `eltype`, we're going to provide a custom row type
Base.eltype(df::cellTable{T}) where {T} = CellRow{T}
Base.length(df::cellTable) = length(getfield(df, :data)[1])
Base.size(df::cellTable) = (length(df),length(getfield(df, :columnNames)))

Base.iterate(df::cellTable, st=0) = st > length(df)-1 ? nothing : (CellRow(st, df), st + 1)

# a custom row type; acts as a "view" into a row of an AbstractVecOrMat
struct CellRow{T} <: AbstractRow
    row::Int
    source::cellTable{T}
end
# required `AbstractRow` interface methods (same as for `AbstractColumns` object before)
# but this time, on our custom row type
getcolumn(r::CellRow, ::Type, col::Int, nm::Symbol) = getfield(getfield(r, :source), :data)[col][getfield(r, :row)]
getcolumn(r::CellRow, i::Int) = getfield(getfield(r, :source), :data)[i][getfield(r, :row)]
getcolumn(r::CellRow, nm::Symbol) = getfield(getfield(r, :source), :data)[getfield(getfield(r, :source), :lookup)[nm]][getfield(r, :row)]
columnnames(df::CellRow) = getfield(getfield(df, :source), :columnNames)

#This is kind of hacky, but it works for now
function Base.show(io::IO, intState::cellTable) 

    data = OffsetArrays.no_offset_view.(getfield(intState,:data))
    df = cellTable(getfield(intState, :columnNames), getfield(intState, :lookup), data)

    pretty_table(df, header_crayon = crayon"yellow bold", vcrop_mode = :middle)
end