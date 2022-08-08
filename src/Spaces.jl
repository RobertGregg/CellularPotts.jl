####################################################
# Space Structure
####################################################

mutable struct CellSpace{N, T<:Integer} <: AbstractSimpleGraph{T}
    ne::T                         #Number of edges
    fadjlist::Vector{Vector{T}}   #Sorted adjacency list [src]: (dst, dst, dst)
    gridSize::NTuple{N, T}        #Size of grid (x,y,z...)
    wrapAround::Bool              #Does the grid wrap around?
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

#3D Version
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

#3D Version
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

function CellSpace(gridSize::NTuple{N, T}; wrapAround=true, cellNeighbors=mooreNeighbors) where {N, T<:Integer}
    
    nodes = prod(gridSize)
    grid = reshape(1:nodes, gridSize...)

    neighborIndices = [cellNeighbors(n) for n in CartesianIndices(grid)]
    fadjlist = [zeros(T, 3^N - 1) for _ in 1:nodes] #assumes worse case senario (i.e. mooreNeighbors)

    #Loop through all neighbor sets and connect the edges
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
        zeros(T,gridSize),
        zeros(T,gridSize))
end


#Allow CellSpace(n,n) in addition to CellSpace((n,n))
CellSpace(gridSize::T...; wrapAround=true, cellNeighbors=mooreNeighbors) where T<:Integer = CellSpace(gridSize; wrapAround, cellNeighbors)

#Needed for induced_subgraph (why?)
function CellSpace{N,T}(n::Integer=0) where {N, T<:Integer}
    fadjlist = [Vector{T}() for _ in one(T):n]
    return CellSpace{N,T}(0, fadjlist, (0,0), true, zeros(T,n,n),zeros(T,n,n))
end

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
    
    wrapType = space.wrapAround ? "Periodic" : "nonPeriodic"
    numNeigbors = maximum(length, space.fadjlist)

    print(io, " $(wrapType) $(numNeigbors)-Neighbor CellSpace{$(N),$(T)}")
end

