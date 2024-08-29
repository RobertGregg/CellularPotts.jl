# # Bringing ODEs to life

# This simulation actually extends [an example](https://diffeq.sciml.ai/latest/features/callback_functions/#Example-3:-Growing-Cell-Population) from the DifferentialEquations.jl documentation describing a growing cell population. Here we take those underlying model dynamics and combine them with CellularPotts.jl

# We begin by loading in both CellularPotts and DifferentialEquations
using CellularPotts, DifferentialEquations

# On the CellularPotts side, we need to create a new CellPotts model which requires a CellSpace, a CellState, and a list of penalties

# The space will use is a 200×200 grid that defaults to periodic boundary conditions
space = CellSpace(200,200)

# In the CellState we specify one epithelial cell with a volume of 200 pixels
# The cell will be positioned at the halfway point within the space. 
state = CellState(
    names = :Epithelial,
    volumes = 200, 
    counts = 1, 
    positions = size(space) .÷ 2
    );

# Keeping the model simple, we'll only include an Adhesion and Volume penalty
penalties = [
    AdhesionPenalty([0 30;
                    30 30]),
    VolumePenalty([5]),
    ]

# Now that we have all the pieces, we can generate a new CPM model.
cpm = CellPotts(space, state, penalties);

# By default, CellularPotts models to not record states as they change overtime to increase computional speed. To have the model record past states we can toggle the appropriate keyword.
cpm.recordHistory = true;

# ## DifferentialEquations.jl setup

# From the DifferentialEquations example, a theoretical protein X was created for each cell that increases linearly in time with rate parameter α
const α = 0.3;

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
    CellDivision!(cpm, cellID)
    return nothing
end

# This will instantiate the ContinuousCallback triggering cell division
ccb = ContinuousCallback(condition,integrator -> affect!(integrator, cpm));

# To pass multiple callbacks into DifferentialEquations, we need to collect them into a set.
callbacks = CallbackSet(pcb, ccb);


# Define the ODE model and solve
u0 = [0.2]
tspan = (0.0,20.0)
prob = ODEProblem(f,u0,tspan)
sol = solve(prob, Tsit5(), callback=callbacks);

# ## Visualization

# We can replicate the plots from the original example   
using Plots, Printf

# Plot the total cell count over time
plot(sol.t, map((x)->length(x),sol[:]), lw=3,
     ylabel="Number of Cells", xlabel="Time", legend=nothing)

# Plot ProteinX dynamics for a specific cell
ts = range(0, stop=20, length=1000)
plot(ts, first.(sol.(ts)), lw=3,
    ylabel="Amount of X in Cell 1", xlabel="Time", legend=nothing)

# Finally, we can create an animation of the CPM to see the cells dividing. I've dropped the first few frames because the first cell takes a while to divide.
proteinXConc = zeros(size(space)...)

anim = @animate for t in Iterators.drop(1:cpm.step.counter,5*timeScale)
    currTime = @sprintf "Time: %.2f" t/timeScale

    cpmt = cpm(t)
    currSol = sol((t+1)/timeScale)

    #Map protein concentrations to space, set :Medium to zero
    for (i,cellID) in enumerate(cpmt.space.nodeIDs)
        proteinXConc[i] = iszero(cellID) ? 0.0 : currSol[cellID]
    end
    
    plt = heatmap(
        proteinXConc,
        c = cgrad([:grey90, :grey, :gold], [0.1, 0.6, 0.9]),
        clims = (0,1),
        title=currTime,
        titlelocation=:left,
        titlefontsize = 32)

    visualize!(plt,cpmt; colorby=:none)

end

gif(anim, "BringingODEsToLife.gif", fps = 30) 