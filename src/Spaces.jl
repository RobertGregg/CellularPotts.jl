####################################################
# Space Structure
####################################################

"""
    CellSpace(gridSize::NTuple{N, T}; isPeriodic=true, neighborhood=:moore)
    CellSpace(gridSize::T...; isPeriodic=true, neighborhood=:moore) where T<:Integer
A concrete type that stores the underlying space cells will occupy.

A `CellSpace()` can be created by supplying a tuple of dimensions or multiple arguments for each dimension, e.g. `CellSpace((3,3,4))` or `CellSpace(3,3,4)`. There are two optional keyword arguements:
 - isPeriodic`::Bool`: Determines if the grid has periodic boundaries
 - neighborhood`::Symbol`: Adjacent cells are either a `:moore` or a `:vonNeumann` neighborhood. Moore neighborhoods include adjacent diagonal positions.
"""
mutable struct CellSpace{N, T<:Integer} <: AbstractSimpleGraph{T}
    ne::T                         #Number of edges
    fadjlist::Vector{Vector{T}}   #Sorted adjacency list [src]: (dst, dst, dst)
    gridSize::NTuple{N, T}        #Size of grid (x,y,z...)
    isPeriodic::Bool              #Does the grid wrap around?
    nodeIDs::Array{T,N}           #Cell's ID for each node
    nodeTypes::Array{T,N}         #Cell's type for each node
end

#This might break things in Graphs.jl, let's see
size(g::CellSpace) = g.gridSize

#CellSpaces are not directed
is_directed(::Type{<:CellSpace}) = false
is_directed(g::CellSpace) = false

#Give compiler information about element types
eltype(g::CellSpace{N,T}) where {N,T} = T
edgetype(g::CellSpace) = SimpleEdge{eltype(g)}

#collect(edges(g::CellSpace)) was giving Vector{Any}
eltype(::Type{Graphs.SimpleGraphs.SimpleEdgeIter{CellSpace{N, T}}}) where {N, T} = Edge{T}

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


#Copied from simplegraph to add edges to the graph
function add_edge!(g::CellSpace{N,T}, e::SimpleEdge{T}) where {N,T}
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


#https://en.wikipedia.org/wiki/Strong_product_of_graphs
#Note: It looks like Graphs.jl will eventually add this function
strong_product(g::G, h::G) where G <: AbstractSimpleGraph = union(cartesian_product(g,h),tensor_product(g,h))

#https://en.wikipedia.org/wiki/King%27s_graph
kingsGrid(dims; periodic=true) = mapfoldl(periodic ? cycle_graph : path_graph, strong_product, reverse(dims))

function CellSpace(gridSize::NTuple{N, T}; isPeriodic=true, neighborhood=:moore) where {N, T<:Integer}

    if neighborhood == :vonNeumann
        g = Graphs.grid(gridSize; periodic=isPeriodic)
    elseif neighborhood == :moore
        g = kingsGrid(gridSize; periodic=isPeriodic)
    else
        error("Unknown cell neighbor option, current options are :moore and :vonNeumann")
    end

    return CellSpace{N,T}(
        ne(g),
        g.fadjlist,
        gridSize,
        isPeriodic,
        zeros(T,gridSize),
        zeros(T,gridSize)
        )
end


####################################################
# Constructors
####################################################

#Allow CellSpace(n,n) in addition to CellSpace((n,n))
CellSpace(gridSize::T...; isPeriodic=true, neighborhood=:moore) where T<:Integer = CellSpace(gridSize; isPeriodic, neighborhood)

#Needed for induced_subgraph (why?)
function CellSpace{N,T}(n::Integer=0) where {N, T<:Integer}
    fadjlist = [Vector{T}() for _ in one(T):n]
    return CellSpace{N,T}(0, fadjlist, (n,n), true, zeros(T,n,n), zeros(T,n,n))
end

####################################################
# Misc functions
####################################################

arrayids(space::CellSpace{N,T}) where {N,T} = space.nodeIDs
arraytypes(space::CellSpace{N,T}) where {N,T} = space.nodeIDs

####################################################
# Show method
####################################################
function show(io::IO, space::CellSpace{N,T}) where {N,T} 
    
    
    for (i,dim) in enumerate(space.gridSize)
        if i < length(space.gridSize)
            print(io, "$(dim)×")
        else
            print(io, "$(dim)")
        end
    end
    
    wrapType = space.isPeriodic ? "Periodic" : "nonPeriodic"
    numNeigbors = maximum(length, space.fadjlist)

    print(io, " $(wrapType) $(numNeigbors)-Neighbor CellSpace{$(N),$(T)}")
end

#https://discourse.julialang.org/t/weird-base-show-behavior-with-custom-struct/91321
show(io::IO, ::MIME"text/plain", space::CellSpace{N,T}) where {N,T}  = show(io, space)