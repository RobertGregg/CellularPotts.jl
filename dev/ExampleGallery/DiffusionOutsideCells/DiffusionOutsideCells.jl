# # Diffusion Outside Cells

# In this example we take a simple approach to allow cells to continuously secrete a species into extracellular space. We assume the concentration of this species is constant within the cell and that it degrades over time after leaving the cell. 

using CellularPotts, DifferentialEquations
using Graphs
using Plots, Printf

# ## Defining the Cellular Potts Model

# We begin by defining the dimensionality of the space cells will occupy (an N×N grid) and a variable to represent the constant concentration value inside each cell.
const N = 200
const p0 = 50.0;

# Next we generate a `CellPotts()` model with 3 cells approximately 1000 pixels in size with adhesion and volume constraints.
cpm = CellPotts(
    CellSpace(N,N), 
    CellTable([:Epithelial],[1000],[3]),
    [AdhesionPenalty([0 30;30 30]), VolumePenalty([5])]
    );

# Record the model.
cpm.record = true


# Because cell space is represented as a graph, we can use the `laplacian_matrix()` function from the Graphs.jl package to estimate the continuous Laplacian operator that describes diffusion (see these two wikipedia articles for more information: [Graph Laplacians](https://en.wikipedia.org/wiki/Discrete_Laplace_operator#Graph_Laplacians), [Laplacian matrix](https://en.wikipedia.org/wiki/Laplacian_matrix#Laplacian_matrix)). This graphical Laplacian matrix can be easily calculated as the difference between the degree matrix and the adjacency matrix. For our purposes, `laplacian_matrix()` generates a sparse $N^2 \times N^2$ matrix which we can multiply by our state vector to simulate diffusion. We also need to specify that our graph is bi-directional (i.e. undirected).
const 𝓛 = laplacian_matrix(cpm.space, dir=:both);

# To couple the Cellular Potts and PDE models, we use a periodic callback that stops the ODE solver and updates the CPM. Additionally, this function updates grid values that cells occupy with the constant intracellular cell concentration (p0).
function cpmUpdate!(integrator, cpm)
    ModelStep!(cpm)

    u = integrator.u

    for i in 1:nv(cpm.space)
        if !iszero(cpm.space.nodeIDs[i])
            u[i] = p0
        end
    end
end;


# This timeScale variable controls how often the callback is triggered. Larger timescales correspond to faster cell movement.
timeScale = 100
cb = PeriodicCallback(integrator -> cpmUpdate!(integrator, cpm), 1/timeScale);

# ## Defining the PDE Model

# Here there are two parameters that control the rate a species degradation (k) and the rate of diffusion (D). Note that the graphical Laplacian corresponds to the **negative** continuous Laplacian which is why we negate this term.
function f(du,u,p,t)
    k, D = p
    du .= -k*u - D*𝓛*u
    return nothing
end

# Initialize the grid so that spaces occupied by cells receive the initial intracellular concentration p0. Also note that because we are using the graphical Laplacian matrix, we need to flatten our state grid into a vector.
idx = findall(!iszero, cpm.space.nodeIDs[:])
u0 = zeros(N^2)
u0[idx] .= p0;

# ## Solving the System
# Simply call the ODE solver with a given time span and values for the parameters. Other ODE solvers will be more efficient but here we demonstrate that even a basic non-stiff solver will suffice. 
tspan = (0.0,20.0)
p = [0.01, 10.0]
prob = ODEProblem(f,u0,tspan,p)
sol = solve(prob, Tsit5(), callback=cb);


# ## Animating the Solution
# As we move forward in time we see the species leave the cell and eventually degrade in the extracellular space.
anim = @animate for t in range(1, cpm.step.stepCounter, step=10)
    
    currTime = @sprintf "Time: %.2f" t/timeScale

    plt = heatmap(
        reshape(sol(t/timeScale), N,N)',
        axis=nothing,
        legend=false,
        framestyle = :box,
        size=(600,600),
        clim = (0,50),
        c = :twilight,
        title=currTime,
        titlefontsize = 36,
        xlims=(0.5, N+0.5),
        ylims=(0.5, N+0.5),)

    cellborders!(plt,cpm(t).space)

    plt
end

gif(anim, "DiffusionOutsideCells.gif", fps = 30)
