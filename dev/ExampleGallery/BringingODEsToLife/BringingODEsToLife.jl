# # Bringing ODEs to life

# This simulation actually extends [an example](https://diffeq.sciml.ai/latest/features/callback_functions/#Example-3:-Growing-Cell-Population) from the DifferentialEquations.jl documentation describing a growing cell population. Here we take those underlying model dynamics and combine them with CellularPotts.jl

# We begin by loading in both CellularPotts and DifferentialEquations
using CellularPotts, DifferentialEquations

# On the CellularPotts side, we need to create a new CellPotts model which requires a CellSpace, a CellTable, and a list of penalties

# The space will use is a 200×200 grid that defaults to periodic boundary conditions
space = CellSpace(200,200)

# In the CellTable we specify one epithelial cell with a volume of 200 pixels
initialCellState = CellTable(
    [:Epithelial],
    [200],
    [1]);

# The cell will be positioned at the halfway point within the space. 
positions = [size(space) .÷ 2]
initialCellState = addcellproperty(initialCellState, :positions, positions)


# From the DifferentialEquations example, a theoretical protein X was created for each cell that increases linearly in time with rate parameter α
const α = 0.3;

# ProteinX needs an initial condition which we set to 0.2. Note that for :Medium (grid locations without a cell) we give a concentration of zero. 
u0 = [0.0, 0.2]
initialCellState = addcellproperty(initialCellState, :ProteinX, u0)

# Keeping the model simple, we'll only include an Adhesion and Volume penalty
penalties = [
    AdhesionPenalty([0 30;
                    30 30]),
    VolumePenalty([5]),
    ]

# Now that we have all the pieces, we can generate a new CPM model.
cpm = CellPotts(space, initialCellState, penalties);

# ## DifferentialEquations.jl setup

# Currently by default CellularPotts models to not record states as they change overtime to increase computional speed. To have the model record past states we can toggle the appropriate keyword.
cpm.record = true;

# As ProteinX evolves over time for each cell, the CPM model also needs to step forward in time to try and minimize its energy. To facilitate this, we can use the callback feature from DifferentialEquations.jl. Here specifically we use the `PeriodicCallback` function which will stop the ODE solve at regular time intervals and run some other function for us (Here it will be the `ModelStep!` function). 

function cpmUpdate!(integrator, cpm)
    ModelStep!(cpm)
end


# This timeScale variable controls how often the callback is triggered. Larger timescales correspond to faster cell movement.
timeScale = 100
pcb = PeriodicCallback(integrator -> cpmUpdate!(integrator, cpm), 1/timeScale);

# The ODE functions are taken directly from the DifferentialEquations example. Each cell is given the following differential equation
# ```math
#   \frac{\mathrm{d} X}{\mathrm{d} t} = \alpha X
# ```
function f(du,u,p,t)
    for i in eachindex(u)
      du[i] = α*u[i]
    end
end

# Also coming from the differential equations example, this callback is tiggered whenever a cell's ProteinX is greater than 1. Basically the cell will divide when when the ProteinX concentration is too large.
condition(u,t,integrator) = 1-maximum(u)

function affect!(integrator,cpm)
    u = integrator.u
    resize!(integrator,length(u)+1)
    cellID = findmax(u)[2]
    Θ = rand()
    u[cellID] = Θ
    u[end] = 1-Θ

    #Adding a call to divide the cells in the CPM
    CellDivision!(cpm, cellID-1)
    return nothing
end

# This will instantiate the ContinuousCallback triggering cell division
ccb = ContinuousCallback(condition,integrator -> affect!(integrator, cpm));

# To pass multiple callback into DifferentialEquations we need to collect them into a set.
callbacks = CallbackSet(pcb, ccb);


# Define the ODE model and solve
tspan = (0.0,20.0)
prob = ODEProblem(f,u0,tspan)
sol = solve(prob, Tsit5(), callback=callbacks);

# ## Visualization

# We can replicate the plots from the original example   
using Plots, Printf, ColorSchemes

# Plot the total cell count over time
plot(sol.t,map((x)->length(x),sol[:]),lw=3,
     ylabel="Number of Cells",xlabel="Time",legend=nothing)

# Plot ProteinX dynamics for a specific cell
ts = range(0, stop=20, length=100)
plot(ts,map((x)->x[2],sol.(ts)),lw=3, ylabel="Amount of X in Cell 1",xlabel="Time",legend=nothing)

# Finally, we can create an animation of the CPM to see the cells dividing. I've dropped the first few frames because the first cell takes a while to divide.
anim = @animate for t in Iterators.drop(1:cpm.step.stepCounter,5*timeScale)
    currTime = @sprintf "Time: %.2f" t/timeScale

    space = cpm(t).space
    
    plotObject = heatmap(
        space.nodeIDs',
        axis=nothing,
        legend = :none,
        framestyle = :box,
        aspect_ratio=:equal,
        size=(600,600),
        c = cgrad(:tol_light, rev=true),
        title=currTime,
        titlefontsize = 36,
        xlims=(0.5, size(space.nodeIDs,1)+0.5),
        ylims=(0.5, size(space.nodeIDs,2)+0.5))

    cellborders!(plotObject,space) #very slow

    plotObject
end

gif(anim, "BringingODEsToLife.gif", fps = 30)
