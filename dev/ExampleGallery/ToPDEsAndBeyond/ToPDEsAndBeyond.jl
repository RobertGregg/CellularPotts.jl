#TODO make callback for CellPotts
#How to deal with the irregular grid?

using CellularPotts, DifferentialEquations
using Graphs
using Plots, Printf

const N = 100
const ΔR = zeros(N,N)
const ΔP = zeros(N,N)
const ΔX = zeros(N,N)

cpm = CellPotts(
    CellSpace(N,N), 
    CellTable([:Epithelial],[500],[1]), 
    [AdhesionPenalty([0 30;30 30]), VolumePenalty([5])]
    );


function cpmUpdate!(integrator, cpm)
    
    #Unlike the "BringingODEsToLife" Model, CPM steps will effect the integrator
    #Here we'll manually perform ModelStep!
    for _ in 1:nv(cpm.space)

        MHStep!(cpm)

        #Check if the cells have moved
        if cpm.step.success
            #This means the target cell got smaller and the source got larger

            #Unpack the states from the integrator
            u = integrator.u
            @views begin
                R = u[:,:,1]
                P = u[:,:,2]
                X = u[:,:,3]
            end
           
            targetCellID = cpm.step.targetCellID
            #Redistribute the lost material into the target cell (skip :Medium)
            if !iszero(targetCellID)
                targetNodes = findall(isequal(targetCellID), cpm.space.nodeIDs)

                totalR = sum(R[targetNodes])
                totalP = sum(P[targetNodes])
                totalX = sum(X[targetNodes])
                
                R[targetNodes] .*= iszero(totalR) ? 1.0 : 1.0 + R[cpm.step.targetNode]/totalR
                P[targetNodes] .*= iszero(totalP) ? 1.0 : 1.0 + P[cpm.step.targetNode]/totalP
                X[targetNodes] .*= iszero(totalX) ? 1.0 : 1.0 + X[cpm.step.targetNode]/totalX
            end

            sourceCellID = cpm.step.sourceCellID
            if !iszero(sourceCellID)
                #Copy over the value from the source to the new cell position
                R[cpm.step.targetNode] = R[cpm.step.sourceNode]
                P[cpm.step.targetNode] = P[cpm.step.sourceNode]
                X[cpm.step.targetNode] = X[cpm.step.sourceNode]

                #Redistribute to account for the added in the source cell 
                sourceNodes = findall(isequal(sourceCellID), cpm.space.nodeIDs)

                totalR = sum(R[sourceNodes])
                totalP = sum(P[sourceNodes])
                totalX = sum(X[sourceNodes])

                R[sourceNodes] .*= iszero(totalR) ? 1.0 : 1.0 - R[cpm.step.targetNode]/totalR
                P[sourceNodes] .*= iszero(totalP) ? 1.0 : 1.0 - P[cpm.step.targetNode]/totalP
                X[sourceNodes] .*= iszero(totalX) ? 1.0 : 1.0 - X[cpm.step.targetNode]/totalX
            else
                #If the source was :Medium, remove the value to conserve mass
                R[cpm.step.targetNode] = 0.0
                P[cpm.step.targetNode] = 0.0
                X[cpm.step.targetNode] = 0.0
            end
        end
    end

    #Increment the step counter
    cpm.step.stepCounter += 1

    return nothing
end
    
while !cpm.step.success
    MHStep!(cpm)
end
    
# This timeScale variable controls how often the callback is triggered. Larger timescales correspond to faster cell movement.
timeScale = 1
cb = PeriodicCallback(integrator -> cpmUpdate!(integrator, cpm), 1/timeScale);


#∇u is updated with Laplacian estimate from u
function ∇²(Δu,u,space)

    #Reset Δu
    Δu .= 0.0

    Δx² = nv(space) #Grid spacing
    D=10.0 #Diffusion coefficient
    h = D/Δx²

    #Loop through vertices skipping any without a cell ID
    for vertex in vertices(space)
        if iszero(space.nodeIDs[vertex])
            continue
        end

        #For each vertex average the neighboring values
        for neighbor in neighbors(space, vertex)
            if space.nodeIDs[vertex] == space.nodeIDs[neighbor]
                @inbounds Δu[vertex] += u[neighbor] - u[vertex]
            end
        end
    end
    
    #Finally scale by grid spacing and diffusion coefficient
    Δu .*= h

    return nothing
end


#Here we begin to initialize the PDE Model
u0 = zeros(N,N,3)

for I in CartesianIndices(cpm.space.nodeIDs)
    if !iszero(cpm.space.nodeIDs[I])
        u0[I,2] = 0.2 + 0.1*rand()
        u0[I,3] = 2.5 + 1.25*rand()
    end
end

tspan = (0.0,200.0)
p = [1.0, 0.1, 1.0, 0.1, 1.0, 0.1, 1.0, 8.0]


function oscillator2D!(du,u,p,t,cpm)

    
    @views begin
        R = u[:,:,1]
        P = u[:,:,2]
        X = u[:,:,3]

        dR = du[:,:,1]
        dP = du[:,:,2]
        dX = du[:,:,3]
    end

    k1,k2,k3,k4,k5,k6,K,n = p

    #Calculate the diffusion
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


prob = ODEProblem( (du,u,p,t) -> oscillator2D!(du,u,p,t,cpm) ,u0,tspan,p)
sol = solve(prob, ROCK4(), callback=cb)



anim = @animate for t in range(tspan...,200)
    currTime = @sprintf "Time: %.2f" t
    heatmap(
        sol(t)[:,:,3],
        axis=nothing,
        framestyle = :box,
        size=(1200,1200),
        clim = (0,5),
        title=currTime,
        titlefontsize = 48)
end

gif(anim, "PDE.gif", fps = 30)


using Statistics

ts = range(tspan...,200)
plot(ts, map((x)->mean(x[:,:,1]),sol.(ts)),
    label="R(t)",
    xlabel = "Time",
    ylabel = "Average Concentration",
    framestyle = :box)
plot!(ts, map((x)->mean(x[:,:,2]),sol.(ts)), label="P(t)")
plot!(ts, map((x)->mean(x[:,:,3]),sol.(ts)), label="X(t)")
