using Revise
using CellularPotts

#Example Model

space = CellSpace(100,100)

initialCellState = newCellState(
    [:Epithelial, :TCells],
    [75, 50],
    [100, 10])

addCellProperty!(initialCellState, :isTumor, false, :TCells)


penalties = [
    AdhesionPenalty([0 20 40;
                    20 90 20;
                    40 20 90]),
    VolumePenalty([5,5])]


cpm = CellPotts(space, initialCellState, penalties);

addCellsRandom!(cpm)



#####################################################################
using GLMakie, Colors

function Edge2Grid(dim)
    gridIndices = LinearIndices(dim)

    x1 = reverse(reshape(gridIndices,dim),dims=1)'[:]
    x2 = circshift(x1,dim[2])

    y1 = reverse(reshape(reverse(gridIndices),dim),dims=2)[:]
    y2 = circshift(y1,dim[1])

    append!(x1,x1[1:dim[1]])
    append!(x2,x2[1:dim[1]])
    append!(y1,y1[1:dim[1]])
    append!(y2,y2[1:dim[1]])

    return [[id1,id2] for (id1,id2) in zip([x1;y1],[x2;y2])]
end



function plotCells(cpm::CellPotts)

    fig = Figure(resolution = (1200, 1200), backgroundcolor = RGBf0(0.98, 0.98, 0.98))
    axSim = fig[1, 1] = Axis(fig, title = "Simulation")

    heatmap!(axSim,
                cpm.visual,
                show_axis = false,
                colormap = :Purples) #:Greys_3
        tightlimits!.(axSim)
        hidedecorations!.(axSim) #removes axis numbers


    edgeConnectors = Edge2Grid(cpm.space.gridSize)
    (m,n) = cpm.space.gridSize

    #Generate all of the edge Connections by putting a point on each cell corner
    horizontal = [Point2f0(x, y) => Point2f0(x+1, y) for x in 0.5:m-0.5, y in 0.5:m+0.5]
    vertical = [Point2f0(x, y) => Point2f0(x, y+1) for y in 0.5:n-0.5, x in 0.5:n+0.5]
    points = vcat(horizontal[:],vertical[:])

    #Determine the transparency of the linesegments
    gridflip = rotl90(cpm.visual) #https://github.com/JuliaPlots/Makie.jl/issues/205

    #Cell borders are outlined in black
    black = RGBA{Float64}(0.0,0.0,0.0,1.0);
    clear = RGBA{Float64}(0.0,0.0,0.0,0.0);

    #Loop through all the grid connected and assign the correct color
    currentEdgeColors = [gridflip[edges[1]]==gridflip[edges[2]] ? clear : black for edges in edgeConnectors];

    linesegments!(
            axSim,
            points,
            color = currentEdgeColors,
            linewidth = 2
        )

    return fig
end


#####################################################################
using Tables
using OffsetArrays

struct cellTable{T <: Vector{OffsetVector{T, AA} where {T, AA<:AbstractVector{T}}}} <: Tables.AbstractColumns
    names::Vector{Symbol}
    lookup::Dict{Symbol, Int}
    data::T
end

# declare that cellTable is a table
Tables.istable(::Type{<:cellTable}) = true
# getter methods to avoid getproperty clash
names(df::cellTable) = getfield(df, :names)
data(df::cellTable) = getfield(df, :data)
lookup(df::cellTable) = getfield(df, :lookup)
# schema is column names and types
Tables.schema(df::cellTable{T}) where {T} = Tables.Schema(names(df), eltype.(data(df)))

# column interface
Tables.columnaccess(::Type{<:cellTable}) = true
Tables.columns(df::cellTable) = df
# required Tables.AbstractColumns object methods
Tables.getcolumn(df::cellTable, ::Type{T}, col::Int, nm::Symbol) where {T} = data(df)[col]
Tables.getcolumn(df::cellTable, nm::Symbol) = data(df)[lookup(df)[nm]]
Tables.getcolumn(df::cellTable, i::Int) = data(df)[i]
Tables.columnnames(df::cellTable) = names(df)


# declare that any cellTable defines its own `Tables.rows` method
Tables.rowaccess(::Type{<:cellTable}) = true
# just return itself, which means cellTable must iterate `Tables.AbstractRow`-compatible objects
Tables.rows(df::cellTable) = df
# the iteration interface, at a minimum, requires `eltype`, `length`, and `iterate`
# for `cellTable` `eltype`, we're going to provide a custom row type
Base.eltype(df::cellTable{T}) where {T} = CellRow{T}
Base.length(df::cellTable) = length(data(df)[1])

Base.iterate(df::cellTable, st=0) = st > length(df)-1 ? nothing : (CellRow(st, df), st + 1)

# a custom row type; acts as a "view" into a row of an AbstractVecOrMat
struct CellRow{T} <: Tables.AbstractRow
    row::Int
    source::cellTable{T}
end
# required `Tables.AbstractRow` interface methods (same as for `Tables.AbstractColumns` object before)
# but this time, on our custom row type
Tables.getcolumn(r::CellRow, ::Type, col::Int, nm::Symbol) = getfield(getfield(r, :source), :data)[col][getfield(r, :row)]
Tables.getcolumn(r::CellRow, i::Int) = getfield(getfield(r, :source), :data)[i][getfield(df, :row)]
Tables.getcolumn(r::CellRow, nm::Symbol) = getfield(getfield(r, :source), :data)[getfield(getfield(r, :source), :lookup)[nm]][getfield(r, :row)]
Tables.columnnames(df::CellRow) = names(getfield(df, :source))



exNames = [:Col1, :Col2, :Col3]
exLookup = Dict(exNames .=> 1:3)

v1 = OffsetVector([1,2,3],0:2)
v2 = OffsetVector([4.0, 5.0, 6.0],0:2)
v3 = OffsetVector(["7","8","9"],0:2)
exData = [v1,v2,v3]

df = cellTable(exNames,exLookup,exData)

using CSV

CSV.write("test.csv",df)