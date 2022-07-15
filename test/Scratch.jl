using Revise
using CellularPotts

#Example Model

space = CellSpace(100,100)

initialCellState = newCellState(
    [:Epithelial, :TCells],
    [75, 50],
    [100, 10])

addCellProperty!(initialCellState, :isTumor, false, :Epithelial)


penalties = [
    AdhesionPenalty([0 20 40;
                    20 90 20;
                    40 20 90]),
    VolumePenalty([5,5])]


cpm = CellPotts(space, initialCellState, penalties);

addCellsRandom!(cpm)

for i=1:10000
    MHStep!(cpm)
end
#####################################################################
# using BenchmarkTools

# @benchmark MHStep!(cpm)

#####################################################################
# using GLMakie, Colors

# function Edge2Grid(dim)
#     gridIndices = LinearIndices(dim)

#     x1 = reverse(reshape(gridIndices,dim),dims=1)'[:]
#     x2 = circshift(x1,dim[2])

#     y1 = reverse(reshape(reverse(gridIndices),dim),dims=2)[:]
#     y2 = circshift(y1,dim[1])

#     append!(x1,x1[1:dim[1]])
#     append!(x2,x2[1:dim[1]])
#     append!(y1,y1[1:dim[1]])
#     append!(y2,y2[1:dim[1]])

#     return [[id1,id2] for (id1,id2) in zip([x1;y1],[x2;y2])]
# end



# function plotCells(cpm::CellPotts)

#     fig = Figure(resolution = (1200, 1200), backgroundcolor = RGBf0(0.98, 0.98, 0.98))
#     axSim = fig[1, 1] = Axis(fig, title = "Simulation")

#     heatmap!(axSim,
#                 cpm.visual,
#                 show_axis = false,
#                 colormap = :Purples) #:Greys_3
#         tightlimits!.(axSim)
#         hidedecorations!.(axSim) #removes axis numbers


#     edgeConnectors = Edge2Grid(cpm.space.gridSize)
#     (m,n) = cpm.space.gridSize

#     #Generate all of the edge Connections by putting a point on each cell corner
#     horizontal = [Point2f0(x, y) => Point2f0(x+1, y) for x in 0.5:m-0.5, y in 0.5:m+0.5]
#     vertical = [Point2f0(x, y) => Point2f0(x, y+1) for y in 0.5:n-0.5, x in 0.5:n+0.5]
#     points = vcat(horizontal[:],vertical[:])

#     #Determine the transparency of the linesegments
#     gridflip = rotl90(cpm.visual) #https://github.com/JuliaPlots/Makie.jl/issues/205

#     #Cell borders are outlined in black
#     black = RGBA{Float64}(0.0,0.0,0.0,1.0);
#     clear = RGBA{Float64}(0.0,0.0,0.0,0.0);

#     #Loop through all the grid connected and assign the correct color
#     currentEdgeColors = [gridflip[edges[1]]==gridflip[edges[2]] ? clear : black for edges in edgeConnectors];

#     linesegments!(
#             axSim,
#             points,
#             color = currentEdgeColors,
#             linewidth = 2
#         )

#     return fig
# end



using TimerOutputs

const to = TimerOutput();

function MHStep!_timed(cpm::CellPotts)

    @timeit to "Search" begin
        #unpack current step structure update
        step = cpm.step


        #Loop through until a good source target is found
        searching = true
        while searching
            #Pick a random location on the graph
            step.sourceNode = rand(1:nv(cpm.space))
            #What cell does it belong to?
            step.sourceCellID = cpm.space.nodeIDs[step.sourceNode]

            #Get all of the unique cell IDs neighboring this Node
            step.neighborNodes = neighbors(cpm.space, step.sourceNode)

            #Choose a target
            step.targetCellID = cpm.space.nodeIDs[rand(step.neighborNodes)]

            #Some checks before attempting a flip
                #In the middle of a cell
                inMiddle = checkMiddle(cpm, step)
                #target is the same as source cell 
                isSource = step.targetCellID == step.sourceCellID

            #if all checks pass, attempt flip
            if !(inMiddle | isSource) 
                searching = false
            end
        end    
    end


    @timeit to "Penalty" begin
        #Calculate the change in energy when source node is modified
        #ΔH =  sum(f(cpm, step) for f in cpm.parameters.penalties)
        ΔH =  applyPenalties(cpm)
    end

    @timeit to "Accept" begin
        #Calculate an acceptance ratio
        acceptRatio = min(1.0,exp(-ΔH/cpm.temperature))


        if rand() < acceptRatio #If the acceptance ratio is large enough
            #Need to update all cell and graph properties
            #---Cell properties---

            #Update cell volumes
            cpm.currentState.volumes[step.sourceCellID] -= 1
            cpm.currentState.volumes[step.targetCellID] += 1

            #TODO update cell perimeters

            #---Graph properties---

            #Cell IDs
            cpm.space.nodeIDs[step.sourceNode] = step.targetCellID
            cpm.space.nodeTypes[step.sourceNode] = cpm.currentState.names[step.targetCellID]
            
            #---Overall properties---
            #Update visual
            cpm.visual[step.sourceNode] = step.targetCellID
            cpm.step.stepCounter += 1
        end
    end

    return nothing
end

function checkMiddle(cpm, step)
    
    for neighbor in step.neighborNodes
        if step.sourceCellID ≠ cpm.space.nodeIDs[neighbor]
            return false
        end
    end

    return true
end

function applyPenalties(cpm)
    ΔH = 0

    for penalty in cpm.penalties
        ΔH += penalty(cpm)
    end

    return ΔH
end


################################

data = [[1,2,3],[:x,:y,:z]]
colNames = [:A,:B]
lookup = Dict(colNames .=> 1:2)
c = cellTable(colNames, lookup, data)