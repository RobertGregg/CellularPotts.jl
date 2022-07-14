
####################################################
# Space Structure
####################################################

mutable struct CellSpace{N, T<:Integer} <: AbstractSimpleGraph{T}
    ne::T                       #Number of edges
    fadjlist::Vector{Vector{T}} #Sorted adjacency list [src]: (dst, dst, dst)
    gridSize::NTuple{N, T}      #Size of grid (x,y,z...)
    wrapAround::Bool            #Does the grid wrap around
    nodeIDs::Vector{T}          #Cell's ID for each node
    nodeTypes::Vector{Symbol}   #Cell's type for each node
end

#CellSpaces are not directed
is_directed(::Type{<:CellSpace}) = false
is_directed(g::CellSpace) = false

#Give compiler information about element types
eltype(g::CellSpace{N,T}) where {N,T} = T
edgetype(g::CellSpace) = SimpleEdge{eltype(g)}

#collect(edges(g::CellSpace)) was giving Vector{Any}
edges(g::CellSpace{N,T}) where {N,T<:Integer} = (SimpleEdge(i,j) for (i,v) in enumerate(g.fadjlist) for j in v if j>i)

####################################################
# Methods to add edges to the space
####################################################

#Copied from simplegraph to check for edges
function has_edge(g::CellSpace{N,T}, s, d) where {N,T}
    verts = vertices(g)
    (s in verts && d in verts) || return false  # edge out of bounds
    @inbounds list_s = g.fadjlist[s]
    @inbounds list_d = g.fadjlist[d]
    if length(list_s) > length(list_d)
        d = s
        list_s = list_d
    end
    return insorted(d, list_s)
end

#Same as above but here the input is an edge type
function has_edge(g::CellSpace{N,T}, e::SimpleEdge{T}) where {N,T}
    s, d = T.(Tuple(e))
    return has_edge(g, s, d)
end

####################################################
# Get node neighbors
####################################################

#The mod1 function always returns a number between 1 and n
#0 gets mapped to n
#n+1 gets mapped to 1
#This is exactly what is needed for periodic boundaries
#Make mod1 work for CartesianIndex (I is the tuple of indices)
Base.mod1(x::CartesianIndex{N}, i::Int...) where N = CartesianIndex(mod1.(x.I,i))

#Get the 4 neighboring nodes around current position
function vonNeumannNeighbors(J::CartesianIndex{2})
    Δx = CartesianIndex(1,0)
    Δy = CartesianIndex(0,1)

    return [     J-Δy,
            J-Δx,     J+Δx,
                 J+Δy]
end

function vonNeumannNeighbors(J::CartesianIndex{3})
    Δx = CartesianIndex(1,0,0)
    Δy = CartesianIndex(0,1,0)
    Δz = CartesianIndex(0,0,1)

    return [             J-Δy,               
            J-Δz,   J-Δx,     J+Δx,   J+Δz,
                         J+Δy,             ]
end

#Get the 8 neighboring nodes around current position
function mooreNeighbors(J::CartesianIndex{2})
    Δx = CartesianIndex(1,0)
    Δy = CartesianIndex(0,1)

    return [J-Δx-Δy,J-Δy,J+Δx-Δy,
            J-Δx,        J+Δx,
            J-Δx+Δy,J+Δy,J+Δx+Δy]
end

function mooreNeighbors(J::CartesianIndex{3})
    Δx = CartesianIndex(1,0,0)
    Δy = CartesianIndex(0,1,0)
    Δz = CartesianIndex(0,0,1)

    return [J-Δx-Δy-Δz,J-Δy-Δz,J+Δx-Δy-Δz,  J-Δx-Δy,J-Δy,J+Δx-Δy,  J-Δx-Δy+Δz,J-Δy+Δz,J+Δx-Δy+Δz,
            J-Δx-Δz,   J-Δz,      J+Δx-Δz,  J-Δx,           J+Δx,  J-Δx+Δz,   J+Δz,      J+Δx+Δz,
            J-Δx+Δy-Δz,J+Δy-Δz,J+Δx+Δy-Δz,  J-Δx+Δy,J+Δy,J+Δx+Δy,  J-Δx+Δy+Δz,J+Δy+Δz,J+Δx+Δy+Δz]
end

####################################################
# Space constructor
####################################################


#Given n nodes, create an empty graph
function CellSpace(gridSize::NTuple{N, T}; wrapAround=true, cellNeighbors=mooreNeighbors) where {N, T<:Integer}
    nodes = prod(gridSize)
    grid = reshape(1:nodes, gridSize...)

    neighborIndices = vec([cellNeighbors(n) for n in CartesianIndices(grid)])
    fadjlist = [zeros(T, 3^N - 1) for _ in 1:nodes]

    for (i,n) in enumerate(neighborIndices)
        if wrapAround
            fadjlist[i] = grid[mod1.(n, gridSize...)]
        else
            for j in eachindex(gridSize)
                filter!(x-> 0 < x.I[j] ≤ gridSize[j], n)
            end
            fadjlist[i] = grid[n]
        end
    end

    #Because of mod1, list might not be sorted
    sort!.(fadjlist)
    
    return CellSpace{N,T}(
        sum(length,fadjlist) ÷ 2, #edges are listed twice (e.g. 1=>2 and 2=>1)
        fadjlist,
        gridSize,
        wrapAround,
        zeros(T,nodes),
        fill(:Medium,nodes))
end

#Allow CellSpace(n,n) in addition to CellSpace((n,n))
CellSpace(gridSize::T...; wrapAround=true, cellNeighbors=mooreNeighbors) where T<:Integer = CellSpace(gridSize; wrapAround, cellNeighbors)