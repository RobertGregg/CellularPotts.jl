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
scene, layout = layoutscene(outer_padding, resolution = (700, 700), backgroundcolor = RGB(0.98, 0.98, 0.98))

#--------Simulation Screen--------
axSim = layout[1, 1] = LAxis(scene, title = "Simulation")

#--------Buttons--------
layout[2,1] = axButtons = GridLayout(tellwidth = false,halign = :left)
buttonsLabels = ["▶","■"]
axButtons[1,1:length(buttonsLabels)] = [LButton(
    scene,
    label = lab,
    width = 70,
    buttoncolor = :grey) for lab in buttonsLabels] 

#Set the buttons to individual variables
playpause, stop = contents(axButtons)

#--------Sliders--------

slideLabels = ["Temperature (β):","Volume\nRestriction (λ)"]
slideRanges = [1:0.01:4, 0.5:0.01:1.5]

sliderSublayout = GridLayout(valign = :top)
layout[1:length(slideLabels), 2] = sliderSublayout

sliders = [labelslider!(scene, 
                          sl,
                          sr;
                          format = x -> "$(x)") for (sl, sr) in zip(slideLabels, slideRanges)]

sliderAxes = layout[1:2,2] = getfield.(sliders, :layout)

sliderSublayout[:v] = sliderAxes

#supertitle = layout[0, 2] = LText(scene, "Parameters")
#--------Update Simulation Screen with Model--------

#Create a new simulation
    #number of cells and individual volumes
    σ=100
    Vd=vcat(0,rand(40:55,σ-1))
testModel = CellPotts(σ=σ,Vd=Vd)

time = Node(0)

heatmap_node = @lift begin
    currentTime = $time
    for t in 1:50_000
        MHStep!(testModel)
    end
    testModel.grid
end

heatmap!(axSim,
         heatmap_node,
         show_axis = false,
         colormap = :Purples) #:Greys_3

tightlimits!.(axSim)
hidedecorations!.(axSim)


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
    axSim,
    points,
    color = lineColors_node,
    linewidth = 2
)


lift(playpause.clicks) do clicks
    if isodd(clicks)
        playpause.label = "𝅛𝅛"
    else
        playpause.label = "▶"
    end
end

#Temperature
lift(sliders[1][:slider].value) do val
    testModel.β = val
end

#display the scene in a new window
scene

runsim = true
stop.clicks[] = 0
while runsim

    #Is the pause button pushed?
    if playpause.label[] == "▶"
        time[] += 1
    end

    #Has the stop button been pushed?
    if  stop.clicks[] == 1
        runsim = false
    end

    sleep(eps())

end
