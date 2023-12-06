# # Travel Time

# This example demonstrates how to record and analyze outcomes from a cellular potts simulation. 

# Here we will track the mean displacement of cell centers with and without active migration. From this data, we can calculate **mean square displacement** as a function of lag time.

# ```math
#   msd(\tau) = <\Delta r(\tau)^2> = <(r(t+\tau) - r(t))^2>
# ```

# Here r(t) represents the current position of the cell and the brackets indicate an average over time. If a cell is moving randomly, this function will be linear. Cells with directional motion will deflect the curve upward. Let's see if we can replicate these theoretic results. 

# ## Create a model containing a cell with and without directional motion.

using CellularPotts, Plots
using Random, Statistics

#Set a random seed for reproducibility 
Random.seed!(314);

#Models have same space and cell initializations
space = CellSpace(200,200)

#Center the moving cell to avoid edges and put stationary cell off into a corner
positions = [(30,30),(100,100)]

#Cells will have same volume
initialCellState = CellState([:MovingCell, :StationaryCell], [500, 500], [1,1])
initialCellState = addcellproperty(initialCellState, :positions, positions)

#Add a migration penalty to one cell encourage cell movement
penalties = [
    AdhesionPenalty([30 30 30;
                    30 20 30
                    30 30 40]),
    VolumePenalty([30, 30]),
    PerimeterPenalty([0, 5]),
    MigrationPenalty(75, [0, 100], size(space)) # 0 means no migration
    ]

#Generate the model
cpm = CellPotts(space, initialCellState, penalties)

#Record model iterations
cpm.recordHistory = true

#Simulate the model
for _ in 1:1000
    ModelStep!(cpm)
end

# # Calculate Average Cell Position

# Given a `CellSpace` and a cell ID, calculate the cell's center.
function meanPosition(space, n)
    
    totalPixels = 0
    avePos = zeros(2) #x and y positions

    #Loop over all points in space
    for i in axes(space,1)
        for j in axes(space,2)
            #Save positions that match cell ID
            if space[i,j] == n
                avePos[1] += i
                avePos[2] += j
                totalPixels += 1
            end
        end
    end

    return avePos ./ totalPixels
end


#Plot trajectories
p1 = plotcpm(cpm)

trajectoryRandom = zeros(2,cpm.step.counter)
trajectoryDirected = zeros(2,cpm.step.counter)

for i in 1:cpm.step.counter
    trajectoryRandom[:,i] .= meanPosition(cpm(i).nodeIDs, 1)
    trajectoryDirected[:,i] .= meanPosition(cpm(i).nodeIDs, 2)
end

plot!(p1, trajectoryRandom[1,:], trajectoryRandom[2,:])
plot!(p1, trajectoryDirected[1,:], trajectoryDirected[2,:])


# # Plot Mean Squared Displacement

MeanSqDis(cpm, τ, id) = mean([sum(abs2, meanPosition(cpm(i+τ).nodeIDs, id) - meanPosition(cpm(i).nodeIDs, id)) for i in 1:τ:cpm.step.counter-τ])

# Here we choose a lag time of 50 steps (arbitrary)

scatter([MeanSqDis(cpm, τ, 1) for τ in 1:50],
    title = "Random Motion",
    xlabel = "Lag Time",
    ylabel = "Mean Squared Displacement",
    legend = false,
    framestyle = :box)

# And plotting the directed motion

scatter([MeanSqDis(cpm, τ, 2) for τ in 1:50],
    title = "Directed Motion",
    xlabel = "Lag Time",
    ylabel = "Mean Squared Displacement",
    legend = false,
    framestyle = :box)

# Here we see that, as a function of lag time, mean squared displacement for random motion does produce a linear relationship whereas directed motion is deflected upward. 