"""
    recordCPM(
        file::String,
        cpm::CellPotts,
        timestamps = 0:300,
        cmap = ColorSchemes.tol_muted;
        framerate=60,
        kwargs...)
Generates an animation of the CPM model.
"""
function recordCPM(
    file::String,
    cpm::CellPotts,
    timestamps = 0:300,
    cmap = ColorSchemes.tol_muted;
    addlegend=false,
    framerate=60,
    kwargs...)

    fig = Figure(resolution = (1200, 1200), backgroundcolor = RGBf(0.98, 0.98, 0.98))
    axSim = fig[1, 1] = Axis(fig, aspect=1)

    timestep = Observable(0) #will increase by one every step

    # heatmap_node is an array that updates when timestep updates
    heatmap_node = @lift begin
        currentTime = $timestep
        ModelStep!(cpm)
        cpm.space.nodeTypes
    end


    distintCellTypes = countcelltypes(cpm) + 1
    cellColors = cgrad(cmap, distintCellTypes, categorical=true, rev=true)

    heatmap!(axSim,
             heatmap_node,
             colormap = cellColors)
        tightlimits!.(axSim)
        hidedecorations!.(axSim) #removes axis numbers

    if addlegend
        labels = String.(unique(cpm.currentState.names))
        elements = [PolyElement(polycolor = cellColors[i]) for i in 1:length(labels)]
        Legend(fig[1,2], elements, labels, labelsize  = 30)
    end

    edgeConnectors = Edge2Grid(cpm.space.gridSize)
    (m,n) = cpm.space.gridSize

    #Generate all of the edge Connections by putting a point on each cell corner
    horizontal = [Point2f(x, y) => Point2f(x+1, y) for x in 0.5:m-0.5, y in 0.5:m+0.5]
    vertical = [Point2f(x, y) => Point2f(x, y+1) for y in 0.5:n-0.5, x in 0.5:n+0.5]
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
    migrationIndex = findfirst(x->typeof(x)==MigrationPenalty, cpm.penalties)
    if !isnothing(migrationIndex)

        heatmap_Gm = @lift begin
            currentTime = $timestep
            Matrix(cpm.penalties[migrationIndex].nodeMemory)
        end
        
        #transparent colors
        colmap = RGBA.(colormap("Purples"),0.5)

        heatmap!(axSim,
                heatmap_Gm,
                colormap = colmap)
    end
    
    
    record(fig, file, timestamps; framerate, kwargs...) do t
        timestep[] += 1
    end


end

fig = Figure(resolution = (1200, 1200), backgroundcolor = RGBf(0.98, 0.98, 0.98))
axSim = fig[1, 1] = Axis(fig, aspect = 1)

cmap = ColorSchemes.tol_muted
cmap = cgrad(cmap, 5, categorical=true, rev=true)

hm = heatmap!(
    axSim,
    rand(1:5,20,20),
    colormap = cmap
)
tightlimits!.(axSim)
hidedecorations!.(axSim)

labels = ["Cell 1","Cell 2","Cell 3","Cell 4","Cell 5"]
elements = [PolyElement(polycolor = cmap[i]) for i in 1:length(labels)]

Legend(fig[1,2], elements, labels,labelsize  = 30)