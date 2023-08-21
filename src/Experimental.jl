#############################################
# Automatic differentiation
#############################################

#Start with simple central difference approximation

# -Fₓ(x) ≈ ∂H/∂σ(x) ⋅ ∂σ(x)/∂x ≈ (1/2h)⋅( H(σ+dₓσ(x)) - H(σ-dₓσ(x)))

#how do you know what direction the neighbor target is facing?

#Ex 3 by 3 gridS
space = reshape(1:9,3,3)
spaceIndex = CartesianIndices(A)

#=
1  4  7
2  5  8
3  6  9
=#

spaceIndex[5] - spaceIndex[8] #CartesianIndex(0, -1)

Tuple(spaceIndex[5] - spaceIndex[8])


#############################################
# Relative penality contributions
#############################################
using CellularPotts, Plots

cpm = CellPotts(
    CellSpace(50, 50),
    CellState(:Epithelial, 300, 1),
    [AdhesionPenalty(fill(30,2,2)), VolumePenalty([5]), PerimeterPenalty([5])]
)

cpm.recordHistory = true

for i=1:100
    ModelStep!(cpm)
end



plot(stack(last(cpm.history.penalty, 100))')

for i in eachindex(cpm.penalties)
    display(histogram(stack(cpm.history.penalty)[i,:], title=cpm.penalties[i]))
end


rows, columns = size(cpm.space)

plt = heatmap(
        cpm.space.nodeTypes',
        c = cgrad(:tol_light, rev=true),
        grid=false,
        axis=nothing,
        legend=:none,
        framestyle=:box,
        aspect_ratio=:equal,
        size = (600,600),
        xlims=(0.5, rows+0.5),
        ylims=(0.5, columns+0.5),
        clim=(0,2)
            )

cellborders!(plt, cpm.space)


recordCPM("OnPatrol.gif", cpm;
    property = :nodeTypes, frameskip=10, c=:RdBu_3)


#############################################
# Parallelization
#############################################


####################################################
# Metropolis–Hasting Step
####################################################

function MHStep_parallel!(cpm::CellPotts)

    #Reset the success flag
    cpm.step.success = false

    source = cpm.step.source
    target = cpm.step.target

    #Pick a random location on the graph
    source.node = rand(1:nv(cpm.space))
    #What cell does it belong to?
    source.id = cpm.space.nodeIDs[source.node]
    source.type = cpm.space.nodeTypes[source.node]

    #Get all cell IDs neighboring this Node
    source.neighbors = neighbors(cpm.space, source.node)

    #Some space locations are forbidden (see TightSpaces example) 
    if isempty(source.neighbors)
        return false
    end

    #Choose a target
    target.node = rand(source.neighbors)
    target.id = cpm.space.nodeIDs[target.node]
    target.type = cpm.space.nodeTypes[target.node]
    
    #Some checks before attempting a flip

    #target and source same cell? 
    if target.id == source.id
        return false
    end

    #sourceCell surrounded by only sourceCells?
    if all(isequal(source.id, cpm.space.nodeIDs[n]) for n in source.neighbors)
        return false
    end

    
    #Get all cell nodes neighboring target node
    target.neighbors = neighbors(cpm.space, target.node)
    
    #Calculate the change in energy when target node is modified
    ΔH =  calculateΔH(cpm)
    
    
    #Calculate an acceptance ratio
    acceptRatio = min(1.0,exp(-ΔH/cpm.temperature))
    
    
    if rand() < acceptRatio #If the acceptance ratio is large enough

        #Moved into accept loop b/c computationally intensive
        if cpm.checkArticulation && isfragmented(cpm)
            return false
        end

        return true
    end

    return false
end


####################################################
# Model Step
####################################################

function ModelStep_parallel!(cpm::CellPotts)

    modelLock = Threads.SpinLock()
    #Repeat MHStep! for each pixel in the model
    Threads.@threads for _ in 1:nv(cpm.space)
        
        if MHStep_parallel!(cpm)
            Threads.lock(modelLock)
            if cpm.recordHistory
                updateHistory!(cpm)
            end

            #Cell IDs
            cpm.space.nodeIDs[target.node] = source.id
            cpm.space.nodeTypes[target.node] = source.type


            #---Cell properties---
            for i in eachindex(cpm.penalties)
                updateMHStep!(cpm, cpm.penalties[i])
            end

            #Finally update the success flag 
            cpm.step.success = true
            Threads.unlock(modelLock)
        end
    end

    #Increment the step counter
    cpm.step.counter += 1

    #Model updates after MCS
    for i in eachindex(cpm.penalties)
        updateModelStep!(cpm, cpm.penalties[i])
    end

    return nothing
end


function f()
    l = Threads.SpinLock()
    x = 0
    Threads.@threads for i in 1:10^7
        Threads.lock(l)
        x += 1  # this block is executed only in one thread
        Threads.unlock(l)
    end
    return x
end