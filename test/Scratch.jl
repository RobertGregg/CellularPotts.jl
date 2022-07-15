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

MHStep!(cpm)

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