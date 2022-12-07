# Logging/Saving simulations results
#This "works" but is very messy

using DataFrames, StatsBase, Tables
using CellularPotts

#Create a CellTable
initialCellState = CellTable(
    [:Epithelial, :TCell],
    [500, 400],
    [7, 5])

#Create a function that generates a dictionary for each column of the CellTable.
# Each column gets a dataframe that logs changes that occur for that column
function makeHistory(initialCellState)
    
    columnDict = Dict{Symbol, DataFrame}()

    for (property, val) in pairs(parent(initialCellState))
        columnDict[property] = DataFrame(time=Int64[], position=Int64[], value=eltype(val)[])
    end

    return columnDict
end


columnDict = makeHistory(initialCellState)


#Randomnly update the CellTable and log the changes

function randUpdate(initialCellState, columnDict)
    

    numRows, numCols = size(initialCellState)

    colNames = Tables.columnnames(initialCellState)

    for i =1:100

        #Get a random column
        randCol = rand(1:numCols-1)
        randRow = rand(1:numRows-1)
        #Get old value
        oldValue = initialCellState[randRow][randCol]

        #Sample that column
        newValue = rand(initialCellState[colNames[randCol]])

        if oldValue ≠ newValue
            #Update the log
            push!(columnDict[colNames[randCol]], [i,randRow,newValue])
        end

    end
end

randUpdate(initialCellState, columnDict)

for (k,v) in pairs(columnDict)
    @show (k,v)
end


function getState(t, columnDict, initialCellState)
    
    out = deepcopy(initialCellState)
    
    #loop through columns in initialCellState
    for (property, df) in pairs(columnDict)
        for row in eachrow(df)
            if row.time > t
                break
            end
            out[property][row.position] = row.value
        end
    end

    return out
end


getState(20, columnDict, initialCellState)


#############################################

using OrdinaryDiffEq

function fish_stock!(ds, s, p, t)
    max_population, h = p
    ds[1] = s[1] * (1 - (s[1] / max_population)) - h
end

stock = 400.0
max_population = 500.0
min_threshold = 60.0


prob = ODEProblem(fish_stock!, [stock], (0.0, Inf), [max_population, 0.0])
integrator = init(prob, Tsit5(); advance_to_tstop = true)



# We step 364 days with this call.
step!(integrator, 30.0, true)


# Only allow fishing if stocks are high enough
integrator.p[2] = integrator.u[1] > min_threshold ? rand(300:500) : 0.0
# Notify the integrator that conditions may be altered
u_modified!(integrator, true)
# Then apply our catch modifier
step!(integrator, 1.0, true)
# Store yearly stock in the model for plotting
stockHistory = integrator.u[1]
# And reset for the next year
integrator.p[2] = 0.0
u_modified!(integrator, true)

#############################################

using OrdinaryDiffEq



function cellCycle!(du, u, p, t)
    α₁, α₂, α₃, β₁, β₂, β₃, K₁, K₂, K₃, n₁, n₂, n₃ = p

    CDK1, PIK1, APC = u

    du[1] = dCDK1 = α₁ - β₁ * CDK1 * APC^n₁ / (K₁^n₁ + APC^n₁)
    du[2] = dPIK1 = α₂*(1-PIK1) * CDK1^n₂ / (K₂^n₂ + CDK1^n₂) - β₂ * PIK1
    du[3] = dAPC = α₃*(1-APC) * PIK1^n₃ / (K₃^n₃ + PIK1^n₃) - β₃ * APC

    return nothing
end

u0 = zeros(3)
p = [0.1, 3.0, 3.0,
     3.0, 1.0, 1.0,
     0.5, 0.5, 0.5,
     8.0, 8.0, 8.0]
t0 = (0.0, 25.0)


prob = ODEProblem(cellCycle!, u0, t0, p)

sol = solve(prob, Tsit5())

#############################################

#= 
Some chemical will reside in a somewhat circular boundary.
Diffusion in and out of the boundary is slow compared to free diffusion
Two reactions for this chemical:
    x --> 2x
    x --> ∅ 
=#
using DifferentialEquations, Graphs, Plots, Printf


dims = (30,30)
numberOfSpecies = 1
numberOfNodes = prod(dims) # number of sites
grid = Graphs.grid(dims)

center = LinearIndices(dims)[dims .÷2...]
starting_state = zeros(Int,numberOfSpecies, numberOfNodes)
starting_state[center] = 25

tspan = (0.0, 10.0)
rates = [2.0, 1.0] # x generation is slower so everthing will be removed eventually

prob = DiscreteProblem(starting_state, tspan, rates)

reactstoch = [[1 => 1],[1 => 1]]
netstoch = [[1 => 1],[1 => -1]]
majumps = MassActionJump(rates, reactstoch, netstoch)

hopConstants = ones(numberOfSpecies, numberOfNodes)
#boundary is harder to hop over
hopConstants[gdistances(grid, center) .== 6] .= -1.0


alg = DirectCRDirect()
jump_prob = JumpProblem(prob,
                        alg,
                        majumps;
                        hopping_constants=hopConstants,
                        spatial_system = grid,
                        save_positions=(true, false))

sol = solve(jump_prob, SSAStepper())


#maxium value obtained
maxu = maximum(maximum.(sol.u))

anim = @animate for t in range(tspan..., 300)
    currTime = @sprintf "Time: %.2f" t
    heatmap(
        reshape(sol(t), dims),
        axis=nothing,
       #clims = (0,maxu),
        framestyle = :box,
        title=currTime)
end

gif(anim, "test.gif", fps = 60)