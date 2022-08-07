
function CellGUI(cpm::CellPotts)

    fig = Figure(resolution = (1200, 1200), backgroundcolor = RGBf(0.98, 0.98, 0.98))
    axSim = fig[1, 1] = Axis(fig, title = "Simulation")

    timestep = Node(1) #will increase by one every step

    # heatmap_node is an array that updates when timestep updates
    heatmap_node = @lift begin
        currentTime = $timestep
        ModelStep!(cpm)
        cpm.space.nodeTypes
    end

    #cmap = ColorSchemes.nipy_spectral
    #cmap = ColorSchemes.linear_bgyw_20_98_c66_n256
    cmap = ColorSchemes.tol_muted
    distintCellTypes = countcelltypes(cpm) + 1

    heatmap!(axSim,
             heatmap_node,
             show_axis = false,
             colormap = cgrad(cmap, distintCellTypes, categorical=true, rev=true))
        tightlimits!.(axSim)
        hidedecorations!.(axSim) #removes axis numbers


    edgeConnectors = Edge2Grid(cpm.space.gridSize)
    (m,n) = cpm.space.gridSize

    #Generate all of the edge Connections by putting a point on each cell corner
    horizontal = [Point2f0(x, y) => Point2f0(x+1, y) for x in 0.5:m-0.5, y in 0.5:m+0.5]
    vertical = [Point2f0(x, y) => Point2f0(x, y+1) for y in 0.5:n-0.5, x in 0.5:n+0.5]
    points = vcat(horizontal[:],vertical[:])

    #Determine the transparency of the linesegments
    gridflip = rotl90(cpm.space.nodeIDs) #https://github.com/JuliaPlots/Makie.jl/issues/205

    #Cell borders are outlined in black
    black = RGBA{Float64}(0.0,0.0,0.0,1.0);
    clear = RGBA{Float64}(0.0,0.0,0.0,0.0);

    #Loop through all the grid connected and assign the correct color
    currentEdgeColors = [gridflip[edges[1]]==gridflip[edges[2]] ? clear : black for edges in edgeConnectors];

    #For each time update, recolor all of the edges
    lineColors_node = @lift begin
        currentTime = $timestep
        
        gridflip = rotl90(cpm.space.nodeIDs)

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

    #Active cell movement
    #TODO use spy()?
    migrationIndex = findfirst(x-> x isa MigrationPenalty, cpm.penalties)
    if !isnothing(migrationIndex)

        heatmap_Gm = @lift begin
            currentTime = $timestep
            Matrix(cpm.penalties[migrationIndex].nodeMemory)
        end
        
        #transparent colors
        colmap = RGBA.(colormap("Purples"),0.5)

        heatmap!(axSim,
                heatmap_Gm,
                show_axis = false,
                colormap = colmap)
    end

    #--------Buttons--------
    #Currently 2 buttons: a play/pause and a stop button
    #Place the buttons below the simulation and align to the left
    fig[2,1] = buttongrid = GridLayout(tellwidth = false, halign = :left)
    buttonsLabels = ["▶","■","Divide!","Kill"]

    #Loop through and assign button labels, width, and color
    buttongrid[1,1:length(buttonsLabels)] = [Button(
        fig,
        label = lab,
        width = 70,
        buttoncolor = :grey) for lab in buttonsLabels] 

    #Set the buttons to individual variables
    playpause, stop, cellDivideButton, cellDeathButton = contents(buttongrid)

    #If the play/pause button is clicked, change the label
    on(playpause.clicks) do clicks
        if isodd(clicks)
            playpause.label = "𝅛𝅛"
        else
            playpause.label = "▶"
        end
    end

    #partition a random cell when button is clicked 
    on(cellDivideButton.clicks) do clicks
        CellDivision!(cpm,rand(1:countcells(cpm)))
    end

    #Choose a random cell to kill
    on(cellDeathButton.clicks) do clicks
        CellDeath!(cpm,rand(1:countcells(cpm)))
    end

    display(fig)

    runsim = true
    stop.clicks[] = 0 #for stop button
    while runsim

        #Is the pause button pushed?
        if playpause.label[] == "▶"
            timestep[] += 1
            notify(timestep)
        end

        #Has the stop button been pushed?
        if  stop.clicks[] == 1
            runsim = false
            GLMakie.destroy!(GLMakie.global_gl_screen()) #close the window
        end

        
        sleep(eps())
    end
end


function Edge2Grid(dim)
    gridIndices = LinearIndices(dim)

    x1 = reverse(reshape(gridIndices,dim),dims=1)'[:]
    x2 = circshift(x1,dim[2])

    y1 = reverse(reshape(reverse(gridIndices),dim),dims=2)[:]
    y2 = circshift(y1,dim[1])

    append!(x1,x1[1:dim[1]])
    append!(x2,x2[1:dim[1]])
    append!(y1,y1[1:dim[1]])
    append!(y2,y2[1:dim[1]])

    return [[id1,id2] for (id1,id2) in zip([x1;y1],[x2;y2])]
end