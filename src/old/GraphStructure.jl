####################################################
# Network Structure
####################################################

#This is mostly just a copy of Graph.jl version of SimpleGraph but now each node was two extra values: cell ID and cell type.

#Possible improvement: if we assume the network is a lattice where each node have 8 neighbors, then we could use an array (n×8) instead of an adjacency list

mutable struct network{T<:Integer} <: AbstractSimpleGraph{T}
    ne::T #number of edges
    fadjlist::Vector{Vector{T}} #Sorted adjacency list [src]: (dst, dst, dst)
    nodeIDs::Vector{T} #Cell's ID for each node
    nodeTypes::Vector{Symbol} #Cell's type for each node
end

#Given n nodes, create an empty graph
function network{T}(n::Integer=0) where T <: Integer
    fadjlist = [Vector{T}() for _ in one(T):n]
    return network{T}(0, fadjlist, zeros(Int,n), fill(:Medium,n))
end

#When generating a network given n number of nodes, provide type information
#network(n::T) where T <: Integer = network{T}(n)

#Generate an empty network
#network() = network{Int}()

#Copied from simplegraph to add edges to the graph
function add_edge!(g::network{T}, e::SimpleEdge{T}) where T
    s, d = T.(Tuple(e))
    verts = vertices(g)
    (s in verts && d in verts) || return false  # edge out of bounds
    @inbounds list = g.fadjlist[s]
    index = searchsortedfirst(list, d)
    @inbounds (index <= length(list) && list[index] == d) && return false  # edge already in graph
    insert!(list, index, d)

    g.ne += 1
    s == d && return true  # selfloop

    @inbounds list = g.fadjlist[d]
    index = searchsortedfirst(list, s)
    insert!(list, index, s)
    return true  # edge successfully added
end


#Create a network given an adjacency matrix
network(adjmx::AbstractMatrix) = network{Int}(adjmx, zeros(Int,size(adjmx,1)), fill(:Medium,size(adjmx,1)))

#Create a network given an adjacency matrix, cell IDs, and types
network(adjmx::AbstractMatrix,  nodeIDs::AbstractVector, nodeTypes::AbstractVector) = network{Int}(adjmx, nodeIDs, nodeTypes)

# Typed to allow for abstraction, e.g. Graph{UInt8}(adjmx)
function network{T}(adjmx::AbstractMatrix, nodeIDs::AbstractVector, nodeTypes::AbstractVector) where T <: Integer
    dima, dimb = size(adjmx)
    isequal(dima, dimb) || throw(ArgumentError("Adjacency / distance matrices must be square"))
    issymmetric(adjmx) || throw(ArgumentError("Adjacency / distance matrices must be symmetric"))

    g = network(T(dima))
    @inbounds for i in findall(triu(adjmx) .!= 0)
        add_edge!(g, i[1], i[2])
    end

    g.nodeIDs = nodeIDs
    g.nodeTypes = nodeTypes
    return g
end

#Graphs are not directed
is_directed(::Type{<:network}) = false

#Give compiler information about element types
eltype(g::network{T}) where T = T
edgetype(::network{T}) where T <: Integer = SimpleEdge{T}

#Copied from simplegraph to check for edges
function has_edge(g::network{T}, s, d) where T
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
function has_edge(g::network{T}, e::SimpleEdge{T}) where T
    s, d = T.(Tuple(e))
    return has_edge(g, s, d)
end


####################################################
# Generate grid network
####################################################

#The mod1 function always returns a number between 1 and n
#0 gets mapped to n
#n+1 gets mapped to 1
#This is exactly what is needed for periodic boundaries
#Make mod1 work for CartesianIndex (I is the tuple of indices)
Base.mod1(x::CartesianIndex{N}, i::Int...) where N = CartesianIndex(mod1.(x.I,i))

#Get the 8 neighboring nodes around current position
function mooreNeighbors(J::CartesianIndex{2}, gridSize::NTuple{2, Int})
    Δx = CartesianIndex(1,0)
    Δy = CartesianIndex(0,1)
    
    return mod1.([J-Δx-Δy,J-Δy,J+Δx-Δy,
                  J-Δx,        J+Δx,
                  J-Δx+Δy,J+Δy,J+Δx+Δy], gridSize...)
end

function mooreNeighbors(J::CartesianIndex{3}, gridSize::NTuple{3, Int})
    Δx = CartesianIndex(1,0,0)
    Δy = CartesianIndex(0,1,0)
    Δz = CartesianIndex(0,0,1)
    
    return mod1.([
        J-Δx-Δy-Δz,J-Δy-Δz,J+Δx-Δy-Δz,
        J-Δx-Δz,   J-Δz,   J+Δx-Δz,
        J-Δx+Δy-Δz,J+Δy-Δz,J+Δx+Δy-Δz,
        
        J-Δx-Δy,J-Δy,J+Δx-Δy,
        J-Δx,        J+Δx,
        J-Δx+Δy,J+Δy,J+Δx+Δy,
        
        J-Δx-Δy+Δz,J-Δy+Δz,J+Δx-Δy+Δz,
        J-Δx+Δz,   J+Δz,   J+Δx+Δz,
        J-Δx+Δy+Δz,J+Δy+Δz,J+Δx+Δy+Δz], gridSize...)
end

#Given a network, add in all the needed edges
function connectGraph!(g::network{Int}, gridSize::NTuple{N, Int}) where N

    grid = reshape(1:prod(gridSize), gridSize...)

    for J in CartesianIndices(grid)
        for K in mooreNeighbors(J,gridSize)
            add_edge!(g, grid[J], grid[K])
        end
    end

    return g
end