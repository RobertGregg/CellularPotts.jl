####################################################
# Space Structure
####################################################

"""
    CellSpace(gridSize::NTuple{N, T}; periodic=true, diagonal=false)
    CellSpace(gridSize::T...; periodic=true, diagonal=false) where T<:Integer
A concrete type that stores the underlying space cells will occupy.

A `CellSpace()` can be created by supplying a tuple of dimensions or multiple arguments for each dimension, e.g. `CellSpace((3,3,4))` or `CellSpace(3,3,4)`. There are two optional keyword arguments:
 - periodic `::Bool`: Determines if the grid has periodic boundaries
 - diagonal `::Bool`: Adjacent cells can have vonNeumann (default) or Moore neighborhoods. VonNeumann neighborhoods **do not** include adjacent diagonal positions.
"""
mutable struct CellSpace{T<:Integer, C, N} <: AbstractSimpleGraph{T}
    ne::T                         #Number of edges
    fadjlist::Vector{Vector{T}}   #Sorted adjacency list [src]: [dst, dst, dst]
    gridSize::NTuple{N, T}        #Size of grid (x,y,z...)
    periodic::Bool                #Does the grid wrap around?
    nodeIDs::Array{T,N}           #Cell's ID for each node
    nodeTypes::Array{T,N}         #Cell's type for each node
end

#This might break things in Graphs.jl, let's see
size(g::CellSpace) = g.gridSize

#CellSpaces are not directed
is_directed(::Type{<:CellSpace}) = false
is_directed(g::CellSpace) = false

#Give compiler information about element types
eltype(g::CellSpace{T,C,N}) where {T,C,N} = T
edgetype(g::CellSpace) = SimpleEdge{eltype(g)}

#collect(edges(g::CellSpace)) was giving Vector{Any}
eltype(::Type{SimpleEdgeIter{CellSpace{T,C,N}}}) where {T,C,N} = Edge{T}

#Forward and backwards adjacencies
fadj(g::CellSpace) = g.fadjlist
fadj(g::CellSpace, v::Integer) = g.fadjlist[v]

badj(g::CellSpace) = fadj(g)
badj(g::CellSpace, v::Integer) = fadj(g, v)

#Get the maximum number of interior node neighbors in the graph
maxNeighborCount(space::CellSpace{T,C,N}) where {T,C,N} = C

####################################################
# Methods to add edges to the space
####################################################

#Copied from simplegraph to check for edges
function has_edge(g::CellSpace, s, d)
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
function has_edge(g::CellSpace{T,C,N}, e::SimpleEdge{T}) where {T,C,N}
    s, d = T.(Tuple(e))
    return has_edge(g, s, d)
end


#Copied from simplegraph to add edges to the graph
function add_edge!(g::CellSpace{T,C,N}, e::SimpleEdge{T}) where {T,C,N}
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

####################################################
# Generate graphs using graph products
####################################################
#This honestly feels like witchcraft 

#https://en.wikipedia.org/wiki/Strong_product_of_graphs
#Note: It looks like Graphs.jl will eventually add this function
strong_product(g, h) = union(cartesian_product(g,h), tensor_product(g,h))

#https://en.wikipedia.org/wiki/King%27s_graph
kingsGrid(dims; periodic=true) = mapfoldl(periodic ? cycle_graph : path_graph, strong_product, reverse(dims))

function CellSpace(gridSize::NTuple{N, T}; periodic=true, diagonal=false) where {N, T<:Integer}

    if diagonal
        g = kingsGrid(gridSize; periodic)
    else
        g = Graphs.grid(gridSize; periodic)
    end
    
    #Calculate the maximum number of neighbors an interior node can have
    C = maximum(length, g.fadjlist)

    return CellSpace{T,C,N}(
        ne(g),
        g.fadjlist,
        gridSize,
        periodic ,
        zeros(T,gridSize),
        zeros(T,gridSize)
        )
end

####################################################
# Generate Graph from matrix
####################################################
function CellSpace(locationMatrix::Matrix{T}; periodic=true, diagonal=false) where T<:Integer

    space = CellSpace(size(locationMatrix); periodic, diagonal)

    #For some reason rem_edge! wasn't removing all the edges
    for (i,val) in enumerate(locationMatrix)
        if iszero(val)
            empty!(space.fadjlist[i])
        end
    end
    space.ne = sum(length.(space.fadjlist)) ÷ 2

    return space
end


####################################################
# Additional Constructors
####################################################

#Allow CellSpace(n,n) in addition to CellSpace((n,n))
CellSpace(gridSize::T...; periodic=true, diagonal=false) where T<:Integer = CellSpace(gridSize; periodic, diagonal)

#Needed for induced_subgraph (why?)
function CellSpace{T,C,N}(n::Integer=0) where {T<:Integer,C,N}
    fadjlist = [Vector{T}() for _ in one(T):n]
    return CellSpace{T,C,N}(0, fadjlist, (n,n), true, zeros(T,n,n), zeros(T,n,n))
end

####################################################
# Misc functions
####################################################

nodeIDs(space::CellSpace) = space.nodeIDs
nodeTypes(space::CellSpace) = space.nodeTypes

####################################################
# Show method
####################################################
function show(io::IO, space::CellSpace{T,C,N}) where {T,C,N} 
    
    for (i,dim) in enumerate(space.gridSize)
        if i < length(space.gridSize)
            print(io, "$(dim)×")
        else
            print(io, "$(dim)")
        end
    end
    
    wrapType = space.periodic ? "Periodic" : "nonPeriodic"

    print(io, " $(wrapType) $(maxNeighborCount(space))-Neighbor CellSpace{$(T),$(N)}")
end

#https://discourse.julialang.org/t/weird-base-show-behavior-with-custom-struct/91321
show(io::IO, ::MIME"text/plain", space::CellSpace)  = show(io, space)