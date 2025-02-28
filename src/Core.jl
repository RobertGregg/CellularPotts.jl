####################################################
# Structure for the model
####################################################

"""
    CellPotts(space, initialCellState, penalties)
A data container that holds information to run the cellular potts simulation.

Requires three inputs:
 - `space` -- a region where cells can exist, generated using `CellSpace()`.
 - `state` -- a table where rows are cells and columns are cell properties, generated using `CellState()`.
 - `penalties` -- a vector of penalties to append to the model.
"""
mutable struct CellPotts{T<:Integer, C, N, S<:NamedTuple, U}
    initialSpace::CellSpace{T,C,N}
    space::CellSpace{T,C,N}
    initialState::CellState{S}
    state::CellState{S}
    penalties::Vector{U}
    step::MHStep{T}
    fragment::Articulation{T}
    temperature::Float64
    history::History{T,C,N}
    recordHistory::Bool
    checkArticulation::Bool
end

function CellPotts(space::CellSpace{T,C,N}, state::CellState{S}, penalties::Vector{P}) where {T,C,N,S,P}

    #See https://github.com/JuliaLang/julia/pull/44131 for why Unions are used
    U = Union{typeof.(penalties)...}

    cpm = CellPotts(
        space,
        space,
        state,
        state,
        U[p for p in penalties],
        MHStep(),
        Articulation(nv(space)),
        20.0,
        History(space),
        false,
        true)

    #Position the cells in the model
    if hasproperty(state, :positions)
        positionCells!(cpm)
    else
        positionCellsRandom!(cpm)
    end

    #Now that the cells are added in, reset the initial states/Spaces
    cpm.initialSpace = deepcopy(cpm.space)
    cpm.initialState = deepcopy(cpm.state)

    return cpm
end

####################################################
# Helper functions for CellPotts
####################################################

"""
    countcells(cpm::CellPotts)
Count the number of cells in the model 
"""
countcells(space::CellSpace) = maxmimum(nodeIDs(space))
countcells(df::CellState) = length(df.cellIDs) - 1
countcells(cpm::CellPotts) = countcells(cpm.state)

"""
countcelltypes(cpm::CellPotts)
Count the number of cell types in the model 
"""
countcelltypes(space::CellSpace) = maxmimum(nodeTypes(space))
countcelltypes(df::CellState) = length(unique(df.typeIDs)) - 1
countcelltypes(cpm::CellPotts) = countcelltypes(cpm.state)

nodeIDs(cpm::CellPotts) = nodeIDs(cpm.space)
nodeTypes(cpm::CellPotts) = nodeTypes(cpm.space)

#Given a node, get all neighboring nodes that match the node's id
cellneighbors(space,node) = Iterators.filter(n->space.nodeIDs[n]==space.nodeIDs[node], neighbors(space,node))
    

####################################################
# Override Base.show
####################################################

function show(io::IO, cpm::CellPotts) 
    println(io,"Cell Potts Model:")

    #Grid
    print(io, "Grid: ")
    for (i,dim) in enumerate(cpm.space.gridSize)
        if i < length(cpm.space.gridSize)
            print(io, "$(dim)×")
        else
            println(io, "$(dim)")
        end
    end

    #Cells and types
    cellCounts = encode(cpm.state.names)
    print(io,"Cell Counts:")
    for (key, value) in cellCounts #remove medium
        if key ≠ :Medium
            print(io," [$(key) → $(value)]")
        end
    end

    if length(cellCounts) > 2 #Includes Medium
        println(io," [Total → $(length(cpm.state.names)-1)]")
    else
        print(io,"\n")
    end

    print(io,"Model Penalties:")
    for p in penaltyType.(cpm.penalties)
        print(io," ", p)
    end

    print(io,"\n")
    println(io,"Temperature: ", cpm.temperature)
    print(io,"Steps: ", cpm.step.counter)
end