#= Possible components for the gui
    - window displaying the simulation
    - ⏵/⏸ button for the simulation 
    - slider for variables (temperature, adhesiveness, etc)
    - display for much of MH steps
    - count for number of cells
=#

using CellularPotts
using AbstractPlotting
using AbstractPlotting.MakieLayout #Is this needed?
using Colors

outer_padding = 30
scene, layout = layoutscene(outer_padding, resolution = (1200, 700), backgroundcolor = RGB(0.98, 0.98, 0.98))

ax1 = layout[1, 1] = LAxis(scene, title = "Simulation")
#playpause = layout[1,2] = LButton(scene, label = "▶",width = 30)

#Create a number simulation
testModel = CellPotts()

time = Node(0)

heatmap_node = @lift begin
    currentTime = $time
    for t in 1:10_000
        MHStep!(testModel)
    end
    testModel.grid
end

heatmap!(ax1,
         heatmap_node,
         show_axis = false,
         colormap = :Greys_3)

tightlimits!.(ax1)
hidedecorations!.(ax1)


#Generate all the adjacent squares in the grid
edgeConnectors = Edge2Grid(testModel.n)

#Generate all of the edge Connections
horizontal = [Point2f0(x, y) => Point2f0(x+1, y) for x in 0:testModel.n-1, y in 0:testModel.n]
vertical = [Point2f0(x, y) => Point2f0(x, y+1) for y in 0:testModel.n-1, x in 0:testModel.n]
points = vcat(horizontal[:],vertical[:])

#Determine the transparency of the linesegments
gridflip = rotl90(testModel.grid) #https://github.com/JuliaPlots/Makie.jl/issues/205

black = RGBA{Float64}(0.0,0.0,0.0,1.0);
clear = RGBA{Float64}(0.0,0.0,0.0,0.0);

currentEdgeColors = [gridflip[edges[1]]==gridflip[edges[2]] ? clear : black for edges in edgeConnectors];

lineColors_node = @lift begin
    currentTime = $time
    
    gridflip = rotl90(testModel.grid)

    for (i,edges) in enumerate(edgeConnectors)
        currentEdgeColors[i] = gridflip[edges[1]]==gridflip[edges[2]] ? clear : black
    end

    currentEdgeColors
end

linesegments!(
    ax1,
    points,
    color = lineColors_node,
    linewidth = 2
)

for i=1:1_000
    time[] = i
    sleep(1e-3)
end

#############################################
xx = 0:0.2:4pi
line1 = lines!(ax1, sin.(xx), xx, color = :red)
scat1 = scatter!(ax1, sin.(xx) .+ 0.2 .* randn.(), xx,
    color = (:red, 0.5), markersize = 15px, marker = '■')


ax2 = layout[1, 2] = LAxis(scene, title = "Shifted Cosine")

line2 = lines!(ax2, cos.(xx), pi .+ xx, color = :blue)
scat2 = scatter!(ax2, cos.(xx) .+ 0.2 .* randn.(), pi .+ xx,
    color = (:blue, 0.5), markersize = 15px, marker = '▲')


linkaxes!(ax1, ax2)

hideydecorations!(ax2, grid = false)

ax1.xlabel = "Amplitude"
ax2.xlabel = "Amplitude"
ax1.ylabel = "Time [ms]"

leg = layout[1, end+1] = LLegend(scene,
    [line1, scat1, line2, scat2],
    ["True", "Measured", "True", "Measured"])

layout[2, 1:2] = leg

trim!(layout)
leg.tellheight = true
leg.orientation = :horizontal



hm_axes = layout[1:2, 3] = [LAxis(scene, title = t) for t in ["Low Activity", "High Activity"]]
heatmaps = [heatmap!(ax, i .+ rand(100, 100)) for (i, ax) in enumerate(hm_axes)]

hm_sublayout = GridLayout()
layout[1:2, 3] = hm_sublayout

# there is another shortcut for filling a GridLayout vertically with
# a vector of content
hm_sublayout[:v] = hm_axes

tightlimits!.(hm_axes)
hidedecorations!.(hm_axes)

for hm in heatmaps
    hm.colorrange = (1, 3)
end

cbar = hm_sublayout[:, 2] = LColorbar(scene, heatmaps[1], label = "Activity Level")
cbar.width = 30
cbar.height = Relative(2/3)

supertitle = layout[0, :] = LText(scene, "Plotting with MakieLayout",
    textsize = 30, font = "Noto Sans Bold", color = (:black, 0.25))


    scene, layout = layoutscene(resolution = (1200, 900))

    ax = layout[1, 1] = LAxis(scene)
    
    toggles = [LToggle(scene, active = ac) for ac in [true, false]]
    labels = [LText(scene, lift(x -> x ? "active" : "inactive", t.active))
        for t in toggles]
    
    layout[1, 2] = grid!(hcat(toggles, labels), tellheight = false)