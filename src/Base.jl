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

using DataFrames
using StatsBase

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

df = newCellState([:Macrophage, :TCell], [10,12], [5,8])


function addCellProperty(df::DataFrame, propertyname, defaultValue, cellName)
    
    newProperty = [name == cellName ? defaultValue : missing for name in df.names]

    df[!,propertyname] = newProperty
end