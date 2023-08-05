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
mutable struct CellPotts{T<:Integer, C, N, V<:NamedTuple, U}
    #TODO Create a little initial space/state struct
    initialSpace::CellSpace{T,C,N}
    space::CellSpace{T,C,N}
    initialState::CellState{V}
    state::CellState{V}
    penalties::Vector{U}
    step::MHStep{T}
    fragment::Articulation
    temperature::Float64
    history::History{T,C,N}
    recordHistory::Bool
    checkArticulation::Bool

    function CellPotts(space::CellSpace{T,C,N}, state::CellState{V}, penalties::Vector{P}) where {T,C,N,V,P}

        #See https://github.com/JuliaLang/julia/pull/44131 for why Unions are used
        U = Union{typeof.(penalties)...}

        cpm =  new{N,T,V,U}(
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
        if :positions ∈ keys(states)
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
# Reading and Writing History
####################################################
#TODO History update is complicated still
function updateHistory!(cpm::CellPotts, step::Int, idx::Int, nodeID::Int, nodeType::Int)
    
    push!(cpm.history.step, step)
    push!(cpm.history.idx, idx)
    push!(cpm.history.nodeID, nodeID)
    push!(cpm.history.nodeType, nodeType)

    push!(cpm.history.penalty, zeros(length(cpm.penalties)))
    for i in eachindex(cpm.penalties)
        last(cpm.history.penalty)[i] = addPenalty!(cpm, cpm.penalties[i])
    end

    return nothing
end

updateHistory!(cpm::CellPotts) = updateHistory!(
    cpm,
    cpm.step.counter,
    cpm.step.target.node,
    cpm.step.source.id,
    cpm.state.typeIDs[cpm.step.source.id])


#Given the history, retieve the space at a given time step
function (cpm::CellPotts)(t)
    
    cpm.history.space.nodeIDs .= cpm.initialSpace.nodeIDs
    cpm.history.space.nodeTypes .= cpm.initialSpace.nodeTypes
    
    stepMatches = 1:searchsortedlast(cpm.history.step, t)

    cpm.history.space.nodeIDs[cpm.history.idx[stepMatches]] .= cpm.history.nodeID[stepMatches]
    cpm.history.space.nodeTypes[cpm.history.idx[stepMatches]] .= cpm.history.nodeType[stepMatches]

    return cpm.history.space #could also return nothing
end
    
####################################################
# Helper functions for CellPotts
####################################################

"""
    countcells(cpm::CellPotts)
    countcells(df::CellState)

Count the number of cells in the model 
"""
countcells(cpm::CellPotts) = countcells(cpm.state)

"""
    countcelltypes(cpm::CellPotts)
    countcelltypes(df::CellState)

Count the number of cell types in the model 
"""
countcelltypes(cpm::CellPotts) = countcelltypes(cpm.state)

nodeIDs(cpm::CellPotts) = nodeIDs(cpm.space)
nodeTypes(cpm::CellPotts) = nodeTypes(cpm.space)

#Given a cellID calculate it's perimeter 
function calcuatePerimeter(cpm::CellPotts, cellID::Int)

    perimeter = 0

    #Loop through all of space and count neighbors
    for (node, id) in enumerate(cpm.space.nodeIDs)
        if id ≠ cellID
            continue
        end

        for neighbor in neighbors(cpm.space, node)
            if cpm.space.nodeIDs[neighbor] ≠ cellID
                perimeter += 1
            end
        end
    end

    return perimeter
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
    print(io,"Steps: ", cpm.step.stepCounter)
end