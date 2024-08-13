####################################################
# Reading and Writing History
####################################################
function updateHistory!(cpm)

    counter = cpm.step.counter
    node = cpm.step.target.node
    id = cpm.step.source.id
    type = cpm.step.source.type

    return updateHistory!(cpm,counter,node,id,type)
end

function updateHistory!(cpm,counter,node,id,type)

    push!(cpm.history.counter, counter)
    push!(cpm.history.node, node)
    push!(cpm.history.id, id)
    push!(cpm.history.type, type)

    return nothing
end


#Given the history, reconstruct the space at a given time step
#TODO find a better way to loop through history without restarting every time
#TODO should maybe return a cpm instance?
function (cpm::CellPotts)(t::Integer)

    space = cpm.history.space
    
    #Reset the space history
    space.nodeIDs .= cpm.initialSpace.nodeIDs
    space.nodeTypes .= cpm.initialSpace.nodeTypes
    
    lastMatch = 1:searchsortedlast(cpm.history.counter, t-1)

    space.nodeIDs[cpm.history.node[lastMatch]] .= cpm.history.id[lastMatch]
    space.nodeTypes[cpm.history.node[lastMatch]] .= cpm.history.type[lastMatch]

    return CellPotts(
        cpm.initialSpace,
        space,
        cpm.initialState,
        cpm.state,
        cpm.penalties,
        MHStep(),
        Articulation(nv(space)),
        20.0,
        History(space),
        false,
        true)
end