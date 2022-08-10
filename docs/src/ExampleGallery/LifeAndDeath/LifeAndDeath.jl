# # Life and Death

#This simulation actually extend [an example](https://diffeq.sciml.ai/latest/features/callback_functions/#Example-3:-Growing-Cell-Population) from the DifferentialEquations.jl documentation describing a growing cell population. Here we take those underlying model dynamics and combine them with CellularPotts.jl


using CellularPotts, DifferentialEquations

space = CellSpace(200,200)

initialCellState = CellTable(
    [:Epithelial],
    [100],
    [1]);

# The cell will be positioned at the halfway point within the space. 
positions = [size(space) .÷ 2]
initialCellState = addcellproperty(initialCellState, :positions, positions)


# As per the growing cell population example, we define a theoretical protein X that increases linearly in time with rate parameter α (defined later)

u0 = [0.2]
initialCellState = addcellproperty(initialCellState, :ProteinX, u0)

#Set Medium to zero
initialCellState.ProteinX[0] = 0.0

penalties = [
    AdhesionPenalty([0 30;
                    30 30]),
    VolumePenalty([5]),
    ]


cpm = CellPotts(space, initialCellState, penalties);

positionCells!(cpm)


#Use a periodic callback to increment the CPM model every time step


function cpmUpdate!(integrator, cpm)
    ModelStep!(cpm)
    #Update ProteinX
    #cpm.currentState.ProteinX .= integrator.u
end

pcb = PeriodicCallback(integrator -> cpmUpdate!(integrator, cpm), 1.0)

# That's mostly it for the CellularPotts side of things, now lets adding in the DifferentialEquations portion
const α = 0.3

function f(du,u,p,t)
    for i in eachindex(u)
      du[i] = α*u[i]
    end
end


condition(u,t,integrator) = 1-maximum(u)

function affect!(integrator,cpm)
    u = integrator.u
    resize!(integrator,length(u)+1)
    cellID = findmax(u)[2]
    Θ = rand()
    u[cellID] = Θ
    u[end] = 1-Θ

    #Divide the cells
    CellDivision!(cpm, cellID)

    return nothing
end


ccb = ContinuousCallback(condition,integrator -> affect!(integrator, cpm))

callbacks = CallbackSet(pcb, ccb)


tspan = (0.0,10.0)
prob = ODEProblem(f,u0,tspan)
sol = solve(prob,callback=callbacks)

#Well, it runs 🤷
