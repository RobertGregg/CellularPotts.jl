# # Life and Death

# This simulation actually extends [an example](https://diffeq.sciml.ai/latest/features/callback_functions/#Example-3:-Growing-Cell-Population) from the DifferentialEquations.jl documentation describing a growing cell population. Here we take those underlying model dynamics and combine them with CellularPotts.jl


using CellularPotts, DifferentialEquations

# ## CellularPotts.jl Setup

# Start by defining a space and the characteristics of the cells
space = CellSpace(200,200)

initialCellState = CellTable(
    [:Epithelial],
    [200],
    [1]);

# The single cell will be positioned at the halfway point within the space. 
positions = [size(space) .÷ 2]
initialCellState = addcellproperty(initialCellState, :positions, positions)


# As per the growing cell population example, we define a theoretical protein X that increases linearly in time with rate parameter α
const α = 0.3;

# Include an initial protein concentration for :Medium and set it to zero
u0 = [0.0, 0.2]
initialCellState = addcellproperty(initialCellState, :ProteinX, u0)

# Keeping the model simple, we'll only include an Adhesion and Volume penalty
penalties = [
    AdhesionPenalty([0 30;
                    30 30]),
    VolumePenalty([5]),
    ]


# Declare the model and position the cells in the model
cpm = CellPotts(space, initialCellState, penalties);
positionCells!(cpm)

# ## DifferentialEquations.jl setup

# Use a periodic callback to increment the CPM model every time step

# Currently there isn't a simple method to log states in CellPotts models (work in progress! 🙂). For now, we need to create an extrernal variable to log how the nodeIDs change over time.
spaceLog = [cpm.space.nodeIDs];

# Next we need to step the CPM in time at regular intervals. We can accomplish this during the differential equations solve using a periodic callback. This function will trigger at regular interval according to the time scale. Larger timescales correspond to faster cell movement.

# Step the model and update the space log 
function cpmUpdate!(integrator, cpm, spaceLog)
    ModelStep!(cpm)
    push!(spaceLog, copy(cpm.space.nodeIDs))
end

# Set the timescale and create the callback
timeScale = 100
pcb = PeriodicCallback(integrator -> cpmUpdate!(integrator, cpm, spaceLog), 1/timeScale)

# This is taken directly from the example. Each cell is given the following differential equation
# ```math
#   \frac{\mathrm{d} X}{\mathrm{d} t} = \alpha X
# ```
function f(du,u,p,t)
    for i in eachindex(u)
      du[i] = α*u[i]
    end
end

# Also from the differential equations example. This callback is tiggered whenever a cell's ProteinX is greater than 1
condition(u,t,integrator) = 1-maximum(u)

function affect!(integrator,cpm)
    u = integrator.u
    resize!(integrator,length(u)+1)
    cellID = findmax(u)[2]
    Θ = rand()
    u[cellID] = Θ
    u[end] = 1-Θ

    #Divide the cells in the CPM
    CellDivision!(cpm, cellID-1)
    return nothing
end

ccb = ContinuousCallback(condition,integrator -> affect!(integrator, cpm));


# Collect the callbacks together
callbacks = CallbackSet(pcb, ccb);

# Define the ODE model and solve
tspan = (0.0,20.0)
prob = ODEProblem(f,u0,tspan)
sol = solve(prob, Tsit5(),callback=callbacks);

# We can replicate the plots from the original example
using Plots, Printf

# Plot the total cell count over time
plot(sol.t,map((x)->length(x),sol[:]),lw=3,
     ylabel="Number of Cells",xlabel="Time")

# Plot ProteinX dynamics for a specific cell
ts = range(0, stop=20, length=100)
plot(ts,map((x)->x[2],sol.(ts)),lw=3, ylabel="Amount of X in Cell 1",xlabel="Time")

# Finally, we can create an animation of the CPM to see the cells dividing. I've dropped the first few frames because the first cell takes a while to divide.
anim = @animate for t in Iterators.drop(eachindex(spaceLog),5*timeScale)
    currTime = @sprintf "Time: %.2f" t/timeScale
    heatmap(spaceLog[t], axis=nothing,legend = :none, size=(500,500),title=currTime)
end

gif(anim, "LifeAndDeath.gif", fps = 30)