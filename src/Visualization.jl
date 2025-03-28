"""
    cellborders!(plotObject, space::CellSpace)
Calculate line borders to differentiate cells in a plot.
"""
function cellborders(space::CellSpace{T, C, 2}) where {T,C}

    (rows,columns) = size(space)

    borders = fill((NaN,NaN),3*rows*columns)
    
    borderEnd = 0
    for r in 1:rows
        for c in 1:columns

            #vertical
            if space.nodeIDs[r,c] ≠ space.nodeIDs[r,mod1(c+1,columns)]
                borders[borderEnd+1] = (c+0.5, r-0.5)
                borders[borderEnd+2] = (c+0.5, r+0.5)

                borderEnd += 3
            end

            #horizontal
            if space.nodeIDs[r,c] ≠ space.nodeIDs[mod1(r+1,rows),c]
                borders[borderEnd+1] = (c-0.5, r+0.5)
                borders[borderEnd+2] = (c+0.5, r+0.5)

                borderEnd += 3
            end

        end
    end

    return borders[1:borderEnd]
end



#Small internal function to darken a given color
#Technically this will also lighten colors if factor > 1
function darken(oldcolor::HSLA,factor)
    HSLA(
        oldcolor.h,
        oldcolor.s,
        oldcolor.l * factor,
        oldcolor.alpha)
end
    
#Convert the old color to an HSLA color
darken(oldcolor::Colorant,factor) = darken(convert(HSLA,oldcolor),factor)


#Want to use @userplot Visualize but need to dispatch on 2D and 3D CellPotts
#So here we do what the macro does manually
mutable struct Visualize{T<:CellPotts}
    args
end

"""
    visualize(
        cpm::CellPotts;
        colorby=:type,
        cellcolors=nothing,
        migrationcolors,
        kwargs...)
Generates a visualization of the Cellular Potts Model.

This plotting function builds on Plots.jl, so any possible keywords can be passed on.

The optional `colorby` arguement can be: 
    1. `:type` to color the cells by cell type
    2. `:id` to give each cell a unique color
    3. `:none` to not color the cells 

To change the background (Medium) color, supply a color to `background`

A custom set of colors can be supplied, e.g., `colormap = [:red,:blue,:green]`

`migrationcolors` accepts a colormap when a MigrationPenalty is present (cannot be applied to 3D visualizations currently).
"""
visualize(args...; kw...) = RecipesBase.plot(Visualize{typeof(args[1])}(args); kw...)


"""
    visualize!(
        plt::AbstractPlot
        cpm::CellPotts;
        colorby=:type,
        cellcolors=nothing,
        migrationcolors=nothing,
        kwargs...)
Generates a visualization of the CPM model.

This function will append the given plot with a Cellular Potts Model visualization 
"""
visualize!(args...; kw...) = RecipesBase.plot!(Visualize{typeof(args[1])}(args); kw...)
visualize!(plt::RecipesBase.AbstractPlot, args...; kw...) = RecipesBase.plot!(plt, Visualize{typeof(args[1])}(args); kw...)



#2D plot recipe
@recipe function f(
    h::Visualize{CellPotts{T,C,2,S,U}};
    colorby=:type,
    cellcolors=nothing,
    migrationcolors=range(RGBA(1,1,1,0.0), RGBA(60/255, 4/255, 104/255, 0.8))
    ) where {T,C,S,U}

    #Ensure the arguement is a cellular potts model
    if !(typeof(h.args[1]) <: CellPotts)
        error("Visualize expects the first arguement to be a Cellular Potts Model")
    end

    if colorby ∉ [:type, :id, :none]
        error("colorby should be :type, :id, or :none")
    end

    #Extract the cell potts model
    cpm = h.args[1]

    numTypes = countcelltypes(cpm)
    
    #Get the data we want to color the cells by
    #replacing zeros NaN will put no color in that location
    if colorby == :type
        data = replace(cpm.space.nodeTypes, 0=>NaN)

        if isnothing(cellcolors)
            if numTypes == 1
                cellcolors = cgrad([palette(:tol_light)[1]])
            elseif numTypes ≤ 6
                cellcolors = cgrad(palette(:tol_light)[1:numTypes])
            else
                cellcolors = cgrad(:rainbow1)
            end
        end

    elseif colorby == :id

        data = replace(cpm.space.nodeIDs, 0=>NaN)
        if isnothing(cellcolors)
            cellcolors = cgrad(:rainbow1)
        end

    end
    

    #Default Values 
    size --> (500,500)
    background --> :grey90

    #How big is the space to plot?
    (rows,columns) = size(cpm.space)
    
    #Forced Values
    grid := false
    axis := nothing
    legend := :none
    aspect_ratio := :equal
    framestyle := :box
    xlims := (0.5, columns+0.5)
    ylims := (0.5, rows+0.5)
        

    #Use heatmap to color the cell (if desired)
    if colorby ≠ :none
        @series begin
            seriestype := :heatmap
            background_color_outside := :white
            colormap := cellcolors
            data
        end
    end

    #Plot cell migration if present
    migrationIndex = findfirst(x->x isa MigrationPenalty, cpm.penalties)
    
    if !isnothing(migrationIndex)
        nodeMemory = cpm.penalties[migrationIndex].nodeMemory
        #nodeMemory = replace(nodeMemory, 0=>NaN)

        @series begin
            seriestype := :heatmap
            colormap := migrationcolors
            colorbar_entry := false
            nodeMemory
        end
    end

    #Generate the cell outlines
    borders = cellborders(cpm.space)

    #This work because NaN values lead to missing segments in path
    @series begin
        seriestype := :path
        linecolor := :black
        linewidth := 1
        borders
    end
end

#3D plot recipe
@recipe function f(h::Visualize{CellPotts{T,C,3,S,U}}; colorby=:type, cellcolors=nothing) where {T,C,S,U}

    #Ensure the arguement is a cellular potts model
    if !(typeof(h.args[1]) <: CellPotts)
        error("Visualize expects the first arguement to be a Cellular Potts Model")
    end

    if colorby ∉ [:type, :id, :none]
        error("colorby should be :type, :id, or :none")
    end

    cpm = h.args[1]
    

    if colorby == :type
        data = cpm.space.nodeTypes
        if isnothing(cellcolors)
            cellcolors = [:dodgerblue, :red, :lightgreen, :snow, :orange, :violet] 
        end

    elseif colorby == :id
        data = cpm.space.nodeIDs
        if isnothing(cellcolors)
            colorrange = range(start=0,stop=1,length=countcells(cpm))
            cellcolors = get(colorschemes[:rainbow1], colorrange) #sample colormap
        end

    elseif colorby == :none
        data = @. Int(!iszero(cpm.space.nodeTypes))
        cellcolors = [RGBA(1,1,1,1.0)] #ideally this would be clear, but then all lines disappear
    end

    #Parse color names to color type
    #If cellcolors is already parsed this will do nothing
    cellcolors =  parse.(Colorant,cellcolors)
    
    #Default Values 
    size --> (500,500)
    linecolor --> :black
    proj_type --> :persp

    #Forced Values
    legend := :none
    framestyle := :box
    lims := (0,maximum(size(cpm.space)))

    (xmax,ymax,zmax) = size(cpm.space)

    #8 corners of a cube
    xp = [0, 0, 0, 0, 1, 1, 1, 1] .- 0.5;
    yp = [0, 1, 0, 1, 0, 0, 1, 1] .- 0.5;
    zp = [0, 0, 1, 1, 1, 0, 0, 1] .- 0.5;

    #Loop from furthest corner from view to closest corner
    #Plot only the sides of the cubes we can see
    #Color the faces differently to fake a lighting source
    for x in 1:xmax, y in ymax:-1:1, z in zmax:-1:1

        #Skip areas with not cell
        if !iszero(data[x,y,z])

            #Index to color by
            i = data[x,y,z]

            #Place a cube face if not blocked by adjacent cube
            #Side 
            if iszero(data[mod1(x+1,xmax),y,z])
                @series begin
                    seriestype := :mesh3d
                    connections := [(5,6,7,8)]
                    color := cellcolors[i]
                    xp .+ x, yp .+ y, zp .+ z
                end
            end

            #Front
            if iszero(data[x,mod1(y-1,ymax),z])
                @series begin
                    seriestype := :mesh3d
                    connections := [(1,3,5,6)]
                    color := darken(cellcolors[i],0.4)
                    xp .+ x, yp .+ y, zp .+ z
                end
            end

            #Top
            if iszero(data[x,y,mod1(z+1,zmax)])
                @series begin
                    seriestype := :mesh3d
                    connections := [(3,4,8,5)]
                    color := darken(cellcolors[i],0.8)
                    xp .+ x, yp .+ y, zp .+ z
                end
            end
        end
    end

end

"""
    record(
        cpm::CellPotts;
        file = "output.gif",
        steps = 300,
        framerate = 30,
        skip = 1,
        colorby = :type,
        kwargs...)
Generates an animation of the Cellular Potts Model.

Change the file type by modifying the file name (e.g. `file = "output.mp4`)
"""
function record(cpm::CellPotts; 
    file="output.gif",
    steps=300,
    framerate=30,
    skip=1,
    colorby=:type,
    cellcolors=nothing,
    kwargs...)
    
    anim = @animate for _ in 1:steps
        ModelStep!(cpm)
        visualize(cpm; colorby=colorby, cellcolors=cellcolors, kwargs...)
    end every skip
    
    return gif(anim, file, fps = framerate)
end

#TODO Add a progress bar to record
#TODO Create a playback function