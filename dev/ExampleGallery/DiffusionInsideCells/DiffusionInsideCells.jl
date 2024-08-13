# # Diffusion Inside Cells

# This is a more complex demonstration to showcase one way we can simulate diffusion inside a cell. There are many ways to couple Cellular Potts Models with PDEs, but here we take a simple approach similar to the "BringingODEsToLife.jl" example. The basic strategy is:
# 1. Discretize the diffusion PDE into a system of ODEs.
# 2. At each ODE time-point, calculate how each species diffuses and add that to the ODE dynamics.
# 3. Periodically stop the ODE solver and manually perform a `ModelStep()`.
# 4. If a cell becomes larger/smaller, re-normalize the total amount of species in the cell to conserve mass.
# 5. Continue until the final ODE solver time-point is reached.

# As usual, load in the required packages:
using CellularPotts, DifferentialEquations
using Graphs
using Plots, Printf, Statistics

# At this point you might be wondering "Why do we need Graphs.jl?". The space that cells occupy generated by `CellSpace()` is actually a Graph! This has a number of advantages (e.g., periodic boundary conditions are solved by connecting nodes on opposite boundaries with an edge), but here we need it to loop over neighboring vertices in the graph to recalculate the Laplacian as the cells move and change shape. 

# Here we define the dimensions of the `CellSpace()` and create containers to hold the diffusion calculation for each species.
const N = 100
const ΔR = zeros(N,N)
const ΔP = zeros(N,N)
const ΔX = zeros(N,N);

# ## Creating the CPM Model

# The `CellPotts()` model requires three inputs (space, cell table, and penalties). Here we create an N×N space with one 500 pixel cell that has penalities for adhesion and volume. See the HelloWorld example for more explanation. 
cpm = CellPotts(
    CellSpace(N,N), 
    CellState(:Epithelial, 500, 1, positions = (N,N) .÷ 2),
    [AdhesionPenalty([0 30; 30 30]), VolumePenalty([5])]
    );

# We can tell the cellular potts model to save at each iteration
cpm.recordHistory = true

# This next function looks complicated but has a lot of repeating parts. Essentially this function is called by the ODE solver at regular time intervals to update the CPM model. The remainder of the function ensures that mass is conserved.
function cpmUpdate!(integrator, cpm)
    
    #Unlike the "BringingODEsToLife" Model, CPM steps will effect the integrator
    #Here we'll manually perform ModelStep!() by repeating the Metropolis Hastings step N² times (i.e. the number of grid points or graph vertices). 
    for _ in vertices(cpm.space)

        MHStep!(cpm)

        #Check if the cells have moved
        if cpm.step.success
            #This means the target cell got smaller and the source got larger

            #Unpack the states from the integrator (R=mRNA, P=Protein, X=Inhibitor)
            u = integrator.u
            @views begin
                R = u[:,:,1]
                P = u[:,:,2]
                X = u[:,:,3]
            end
           
            #Extract the target cell ID
            targetCellID = cpm.step.target.id

            #Redistribute the lost material into the target cell (skip if target is not a cell)
            if !iszero(targetCellID)

                #Find all the current locations for the target cell
                targetNodes = findall(isequal(targetCellID), cpm.space.nodeIDs)

                #Calculate the current total mass for each species (lower than it was originally)
                totalR = sum(R[targetNodes])
                totalP = sum(P[targetNodes])
                totalX = sum(X[targetNodes])

                #Ignore if the total is zero, otherwise scale each value in the cell up to account for the loss 
                R[targetNodes] .*= iszero(totalR) ? 1.0 : 1.0 + R[cpm.step.target.node]/totalR
                P[targetNodes] .*= iszero(totalP) ? 1.0 : 1.0 + P[cpm.step.target.node]/totalP
                X[targetNodes] .*= iszero(totalX) ? 1.0 : 1.0 + X[cpm.step.target.node]/totalX
            end

            #Extract the source cell ID
            sourceCellID = cpm.step.source.id

            #Redistribute the gained material into the source cell (skip if source is not a cell)
            if !iszero(sourceCellID)

                #Copy over the value from the source to avoid artificial diffusion gradients 
                R[cpm.step.target.node] = R[cpm.step.source.node]
                P[cpm.step.target.node] = P[cpm.step.source.node]
                X[cpm.step.target.node] = X[cpm.step.source.node]

                #Find all the current locations for the source cell
                sourceNodes = findall(isequal(sourceCellID), cpm.space.nodeIDs)
                
                #Calculate the current total mass for each species (higher than it was originally)
                totalR = sum(R[sourceNodes])
                totalP = sum(P[sourceNodes])
                totalX = sum(X[sourceNodes])
                
                #Redistribute to account for the added in the source cell 
                R[sourceNodes] .*= iszero(totalR) ? 1.0 : 1.0 - R[cpm.step.target.node]/totalR
                P[sourceNodes] .*= iszero(totalP) ? 1.0 : 1.0 - P[cpm.step.target.node]/totalP
                X[sourceNodes] .*= iszero(totalX) ? 1.0 : 1.0 - X[cpm.step.target.node]/totalX

            else
                #If the source was not a cell, remove the value to conserve mass
                R[cpm.step.target.node] = 0.0
                P[cpm.step.target.node] = 0.0
                X[cpm.step.target.node] = 0.0
            end
        end
    end

    #Increment the step counter
    cpm.step.counter += 1

    return nothing
end;

    
# This timeScale variable controls how often the callback is triggered. Larger timescales correspond to faster cell movement.
timeScale = 1;

# Finally we can put the cpm updater into the ODE solver callback
cb = PeriodicCallback(integrator -> cpmUpdate!(integrator, cpm), 1/timeScale);

# ## Discretize the Laplacian for Diffusion

# This function updates Δu with a Laplacian estimate from the current species u.
function ∇²(Δu,u,space)

    #Reset Δu
    Δu .= 0.0

    Δx² = nv(space) #Grid spacing
    D=10.0          #Diffusion coefficient
    h = D/Δx²       #Factor to multiply Δu by

    #Loop through vertices skipping any not apart of a cell
    for vertex in vertices(space)

        #If the nodeID is zero, then no cell is present
        if iszero(space.nodeIDs[vertex])
            continue
        end

        #For each cell vertex average the neighboring values if they are in the cell
        for neighbor in neighbors(space, vertex)
            if space.nodeIDs[vertex] == space.nodeIDs[neighbor]
                @inbounds Δu[vertex] += u[neighbor] - u[vertex]
            end
        end
    end
    
    #Finally scale by the grid spacing and diffusion coefficient
    Δu .*= h

    return nothing
end;


# ## Initialize the PDE Model

# The model we're using is called the "Goodwin Model" which was introduced by B. Goodwin in 1965 [^1]. I discovered this model from a [lecture series](http://mcb111.org/w11/w11-lecture.html) by Dr. Elena Rivas when I was looking for a simple cell signaling model with interesting properties (in this case oscillations). She summarizes the model as follows:

# ```@raw html
#   <p style="text-align:center;">
#   <img src="https://github.com/RobertGregg/CellularPotts.jl/blob/master/docs/src/ExampleGallery/DiffusionInsideCells/goodwin_model.png?raw=true" width="445">
#   </p>
# ```

# ```math
# \begin{aligned}
#  \frac{dR}{dt} = k_1 \frac{K^n}{K^n + X^n} - k_2 R  \\
#  \frac{dP}{dt} = k_3 R - k_4 P \\
#  \frac{dX}{dt} = k_5 P - k_6 X
# \end{aligned}
# ```

# The model consists of three states: an RNA species (R) that produces a Protein (P) which then results in the production of an inhibitor (X). The inhibitor slows down the production of RNA creating a negative feedback loop. 

# This system produces dampened oscillations at values of n<8 and periodic oscillation in species concentrations for n>8.

# ```@raw html
#   <p style="text-align:center;">
#   <img src="https://github.com/RobertGregg/CellularPotts.jl/blob/master/docs/src/ExampleGallery/DiffusionInsideCells/goodwin_oscillations.png?raw=true" width="445">
#   </p>
# ```

# Note that this model isn't very realistic because n is usually interpreted as the degree of cooperatively between the inhibitor and the promoter region (meaning more than 8 molecules would have to bind simultaneously). 

# The only modification we make to this model is to add a diffusion term to each species (ΔR, ΔP, and ΔX). 

# Here we set up the initial condition (3 N×N grids for each state). 
u0 = zeros(N,N,3);

# Cells begin with no RNA, a value of 0.2 for P, and 2.5 for X (all with some noise added)
for I in CartesianIndices(cpm.space.nodeIDs)
    if !iszero(cpm.space.nodeIDs[I])
        u0[I,2] = 0.2 + 0.1*rand()
        u0[I,3] = 2.5 + 1.25*rand()
    end
end

# Set the time span to simulate the model
tspan = (0.0,200.0);

# Give values to all the parameters in the model (k₁, k₂, k₃, k₄, k₅, k₆, K, n)
p = [1.0, 0.1, 1.0, 0.1, 1.0, 0.1, 1.0, 8.0];


# Define the main function to simulate the ODE model
function oscillator2D!(du,u,p,t,cpm)

    #Unpack the model variables
    @views begin
        R = u[:,:,1]
        P = u[:,:,2]
        X = u[:,:,3]

        dR = du[:,:,1]
        dP = du[:,:,2]
        dX = du[:,:,3]
    end

    #Unpack the parameters
    k1,k2,k3,k4,k5,k6,K,n = p

    #Calculate the how the states diffuse in one time-step
    ∇²(ΔR,R, cpm.space)
    ∇²(ΔP,P, cpm.space)
    ∇²(ΔX,X, cpm.space)

    #Goodwin Model R=mRNA, P=Protein, X=Inhibitor
    @. begin
        dR = cpm.space.nodeTypes * k1 * K^n/(K^n+X^n) - k2*R + ΔR
        dP = k3*R - k4*P + ΔP
        dX = k5*P - k6*X + ΔX
    end

    return nothing
end

# Create an ODE model
prob = ODEProblem( (du,u,p,t) -> oscillator2D!(du,u,p,t,cpm), u0, tspan, p);

# Give the callback to the solver and simulate the model
sol = solve(prob, ROCK4(), callback=cb);

# ## Exploring the Solution

# Similar to the "BringingODEsToLife" example, we can animate the PDE solution
anim = @animate for t in range(tspan...,200)
    currTime = @sprintf "Inhibitor X Concentration\nTime: %.2f" t
    heatmap(
        sol(t)[:,:,3],
        colormap=:dense,
        clim = (0,5),
        title=currTime,
        titlelocation=:left,
        background_color_outside = :white,
        titlefontsize = 16)
    visualize!(cpm(floor(Int,t)); colorby=:none)
end

gif(anim, "DiffusionInsideCells.gif", fps = 30)

# Note that you can see the oscillations in the inhibitor over time by the color change.

# Finally we can replicate the figure Dr. Elena Rivas created by averaging the cells concentrations over time

# Evaluate the solution at 200 time-points
ts = range(tspan...,200);

# Plot the mean value of each state
plot(ts, map((x)->mean(x[:,:,1]),sol.(ts)),
    label="R(t)",
    xlabel = "Time",
    ylabel = "Average Concentration",
    framestyle = :box);
plot!(ts, map((x)->mean(x[:,:,2]),sol.(ts)), label="P(t)");
plot!(ts, map((x)->mean(x[:,:,3]),sol.(ts)), label="X(t)")


# [^1]: Goodwin, Brian C. "Oscillatory behavior in enzymatic control processes." Advances in enzyme regulation 3 (1965): 425-437.