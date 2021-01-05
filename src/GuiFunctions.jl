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


function CellGUI(n,cellNumber,cellSize)

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
    slideRanges = [1:0.01:4, 0:10]

    Sublayout = GridLayout(valign = :top,tellheight=false)
    layout[1, 2] = Sublayout

    sliders = [labelslider!(scene, 
                            sl,
                            sr;
                            format = x -> "$(x)") for (sl, sr) in zip(slideLabels, slideRanges)]

    sliderAxes = layout[1:2,2] = getfield.(sliders, :layout)

    Sublayout[:v] = sliderAxes

    #supertitle = layout[0, 2] = LText(scene, "Parameters")
    #--------Update Simulation Screen with Model--------

    #Create a new simulation
    #number of cells and individual volumes
    σ=cellNumber + 1
    Vd=vcat(0,cellSize)
    testModel = CellPotts(n=n,σ=σ,Vd=Vd)

    timestep = Node(1)
    frameskip = 50_000

    heatmap_node = @lift begin
        currentTime = $timestep
        for t in 1:frameskip
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
        currentTime = $timestep
        
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

    #Lagrange multiplier for volume
    lift(sliders[2][:slider].value) do val
        testModel.λ[2:end] .= val #skip the medium
    end

    display(scene)

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
        println(testModel.H)

    end
end


CellGUI(100,50,fill(40,50))
