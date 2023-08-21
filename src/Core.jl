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

    function CellPotts(space::CellSpace{T,C,N}, state::CellState{S}, penalties::Vector{P}) where {T,C,N,S,P}

        #See https://github.com/JuliaLang/julia/pull/44131 for why Unions are used
        U = Union{typeof.(penalties)...}

        cpm =  new{T,C,N,S,U}(
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
        #TODO need a better way to position cells
        if :positions ∈ keys(state)
            positionCells!(cpm)
        else
            positionCellsRandom!(cpm)
        end

        #Now that the cells are added in, reset the initial states/Spaces
        cpm.initialSpace = deepcopy(cpm.space)
        cpm.initialState = deepcopy(cpm.state)

        return cpm
    end
end

####################################################
# Helper functions for CellPotts
####################################################

"""
    countcells(cpm::CellPotts)
Count the number of cells in the model 
"""
countcells(cpm::CellPotts) = countcells(cpm.state)

"""
    countcelltypes(cpm::CellPotts)
Count the number of cell types in the model 
"""
countcelltypes(cpm::CellPotts) = countcelltypes(cpm.state)

nodeIDs(cpm::CellPotts) = nodeIDs(cpm.space)
nodeTypes(cpm::CellPotts) = nodeTypes(cpm.space)

#Given a node, get all neighboring nodes that match the node's id
cellneighbors(space,node) = Iterators.filter(n->space.nodeIDs[n]==space.nodeIDs[node], neighbors(space,node))

####################################################
# Reading and Writing History
####################################################
function updateHistory!(cpm, counter, node, id, type)

    push!(cpm.history.counter, counter)
    push!(cpm.history.position, node)
    push!(cpm.history.nodeID, id)
    push!(cpm.history.nodeType, type)

    return nothing
end

updateHistory!(cpm::CellPotts) = updateHistory!(
    cpm, 
    cpm.step.counter,
    cpm.step.target.node,
    cpm.step.source.id,
    cpm.step.source.type)


#Given the history, reconstruct the space at a given time step
#TODO find a better way to loop through history with restarting every time
function (cpm::CellPotts)(t::Integer)

    space = cpm.history.space
    
    #Reset the space history
    space.nodeIDs .= cpm.initialSpace.nodeIDs
    space.nodeTypes .= cpm.initialSpace.nodeTypes
    
    lastMatch = 1:searchsortedlast(cpm.history.counter, t)

    space.nodeIDs[cpm.history.position[lastMatch]] .= cpm.history.nodeID[lastMatch]
    space.nodeTypes[cpm.history.position[lastMatch]] .= cpm.history.nodeType[lastMatch]

    return space #could also return nothing
end
    

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
    cellCounts = countmap(cpm.state.names)
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