#This is mostly just a copy of Graph.jl version of SimpleGraph but now each node was two extra values: cell ID and cell type.


####################################################
# Network Structure
####################################################

#Possible improvement: if we assume the network is a lattice where each node have 4 neighbors, then we could use an array (n×4) instead of an adjacency list

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
network(n::T) where T <: Integer = network{T}(n)

#Generate an empty network
network() = network{Int}()

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