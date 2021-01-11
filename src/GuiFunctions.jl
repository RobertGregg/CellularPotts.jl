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


function CellGUI(CPM)

    #Create a blank scene with some padding
    outer_padding = 30
    scene, layout = layoutscene(outer_padding, resolution = (700, 700), backgroundcolor = RGB(0.98, 0.98, 0.98))

    #--------Simulation Screen--------
    axSim = layout[1, 1] = LAxis(scene, title = "Simulation")

    #--------Buttons--------
    #Currently 2 buttons: a play/pause and a stop button
    #Place the buttons below the simulation and align to the left
    layout[2,1] = axButtons = GridLayout(tellwidth = false,halign = :left)
    buttonsLabels = ["▶","■"]
    #Loop through and assign button labels, width, and color
    axButtons[1,1:length(buttonsLabels)] = [LButton(
        scene,
        label = lab,
        width = 70,
        buttoncolor = :grey) for lab in buttonsLabels] 

    #Set the buttons to individual variables
    playpause, stop = contents(axButtons)

    #--------Sliders--------
    #The first slider controls the inverse temperature (higher temperatures ⟹ accepting higher energy proposals)
    #The second slider changes volume constraint (lower values ⟹ cells more easily drift from desired volume)
    slideLabels = ["Temperature (β):","Volume\nRestriction (λ)"]
    slideRanges = [1:0.01:4, 1:10]

    #Put the sliders to the right of the simulation (tellheight=false means the slider and simulation heights can differ)
    Sublayout = GridLayout(valign = :top,tellheight=false)
    layout[1, 2] = Sublayout

    #Loop through sliders and assign labels and ranges
    sliders = [labelslider!(scene, 
                            sl,
                            sr;
                            format = x -> "$(x)") for (sl, sr) in zip(slideLabels, slideRanges)]

    #Vertically stack the sliders
    sliderAxes = layout[1:2,2] = getfield.(sliders, :layout)
    Sublayout[:v] = sliderAxes

    #--------Update Simulation Screen with Model--------
    timestep = Node(1) #A Node is a variable being observed by the simulation
    frameskip = 1 #The MHStep will often fail b/c of bad random choice to change grid, so nothing will update

    #Everytime time is updated, run MHStep until the next frameskip
    heatmap_node = @lift begin
        currentTime = $timestep
        for t in 1:frameskip
            MHStep!(CPM)
        end
        CPM.grid
    end

    #Create a heatmap on the simulation axis (layout 1,1)
    heatmap!(axSim,
            heatmap_node,
            show_axis = false,
            colormap = :Purples) #:Greys_3
    tightlimits!.(axSim)
    hidedecorations!.(axSim) #removes axis numbers


    #Generate all the adjacent squares in the grid
    edgeConnectors = Edge2Grid(CPM.n)

    #Generate all of the edge Connections by putting a point on each cell corner
    horizontal = [Point2f0(x, y) => Point2f0(x+1, y) for x in 0:CPM.n-1, y in 0:CPM.n]
    vertical = [Point2f0(x, y) => Point2f0(x, y+1) for y in 0:CPM.n-1, x in 0:CPM.n]
    points = vcat(horizontal[:],vertical[:])

    #Determine the transparency of the linesegments
    gridflip = rotl90(CPM.grid) #https://github.com/JuliaPlots/Makie.jl/issues/205

    #Cell borders are outlined in black
    black = RGBA{Float64}(0.0,0.0,0.0,1.0);
    clear = RGBA{Float64}(0.0,0.0,0.0,0.0);
    
    #Loop through all the grid connected and assign the correct color
    currentEdgeColors = [gridflip[edges[1]]==gridflip[edges[2]] ? clear : black for edges in edgeConnectors];

    #For each time update, recolor all of the edges
    lineColors_node = @lift begin
        currentTime = $timestep
        
        gridflip = rotl90(CPM.grid)

        for (i,edges) in enumerate(edgeConnectors)
            currentEdgeColors[i] = gridflip[edges[1]]==gridflip[edges[2]] ? clear : black
        end

        currentEdgeColors
    end

    #Plot all of the line segments onto the simulation
    linesegments!(
        axSim,
        points,
        color = lineColors_node,
        linewidth = 2
    )

    #If the play/pause button is clicked, change the label
    lift(playpause.clicks) do clicks
        if isodd(clicks)
            playpause.label = "𝅛𝅛"
        else
            playpause.label = "▶"
        end
    end

    #Temperature slider update
    lift(sliders[1][:slider].value) do val
        CPM.β = val
    end

    #Volume constraint slider update
    lift(sliders[2][:slider].value) do val
        CPM.λ[1:end] .= val #skip the medium (with index zero)
    end

    #Active cell movement (patrol)
    if CPM.Ma ≠ 0
        #do I need this?
        heatmap_Gm = @lift begin
            currentTime = $timestep
            CPM.Gm
        end
        
        #transparent colors
        colmap = RGBA.(colormap("Reds"),0.5)

        #Create a heatmap on the simulation axis (layout 1,1)
        heatmap!(axSim,
                heatmap_Gm,
                show_axis = false,
                colormap = colmap) #:Greys_3
        tightlimits!.(axSim)
        hidedecorations!.(axSim) #removes axis numbers
    end    

    #Create a window
    display(scene)

    #Keep running the simulation until stopped
    runsim = true
    stop.clicks[] = 0
    while runsim

        #Is the pause button pushed?
        if playpause.label[] == "▶"
            timestep[] += 1
        end

        #Has the stop button been pushed?
        
        if  stop.clicks[] == 1
            runsim = false
            #GLMakie.destroy!(GLMakie.global_gl_screen()) #this might work to close window when stop is pressed?
        end

        sleep(eps())

    end
end


#Create a new simulation
CPM = CellPotts(n=50,σ=1,Vd=[400])

for i=1:1_000_000
    MHStep!(CPM)
    i % 10000 == 0 ? println(i) : nothing
end


CellGUI(CPM)
