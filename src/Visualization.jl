"""
    cellborders!(plotObject, space::CellSpace)
Add line borders to differentiate cells in a plot.
"""
function cellborders!(plotObject, space::CellSpace{2, T}) where T

    (row,col) = size(space)


    x = fill(NaN, 3*row*col)
    y = fill(NaN, 3*row*col)
    
    count = 0
    for r in 1:row
        for c in 1:col

            #vertical
            if space.nodeIDs[r,c] ≠ space.nodeIDs[r,mod1(c+1,col)]
                x[count+1] = r-0.5
                x[count+2] = r+0.5

                y[count+1] = c+0.5
                y[count+2] = c+0.5
                count += 3
            end

            #horizontal
            if space.nodeIDs[r,c] ≠ space.nodeIDs[mod1(r+1,row),c]
                x[count+1] = r+0.5
                x[count+2] = r+0.5

                y[count+1] = c-0.5
                y[count+2] = c+0.5
                count += 3
            end

        end
    end

    plot!(plotObject, x[1:count], y[1:count],color=:black, lw=1, legend=false)

    return nothing
end

#TODO Need proper scale
cellMovement!(plotObject, cpm) = cellMovement!(plotObject, cpm, 1)

function cellMovement!(plotObject, cpm, colorMax)

    #Active cell movement
    migrationIndex = findfirst(x->x isa MigrationPenalty, cpm.penalties)

    if isnothing(migrationIndex)
        return nothing
    end
        
    #transparent colors
    transparentToColor = range(colorant"rgba(255,255,255,0.0)", colorant"rgba(60,4,104,0.6)")

    normConst = colorMax / cpm.penalties[migrationIndex].maxAct

    heatmap!(plotObject,
            normConst .* Matrix(cpm.penalties[migrationIndex].nodeMemory)',
            c = transparentToColor,
            colorbar_entry = false)

    return nothing
end

#TODO Color issue with Patrol, set lower clim (?)
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
    cpm::CellPotts{2, T, V, U};
    timestamps = 0:3000,
    c = cgrad(:tol_light, rev=true),
    figureSize = (600,600),
    property = :nodeIDs,
    legend=:none,
    framerate=30,
    frameskip = 1,
    kwargs...) where {T,V,U}

    (rows,columns) = size(cpm.space)

    #TODO How to color by any property?
    if property == :nodeIDs
        colorMax = countcells(cpm)
    elseif property == :nodeTypes
        colorMax = countcelltypes(cpm)
    else
        colorMax = Inf
    end

    anim = @animate for t in timestamps

        plotObject = heatmap(
            getproperty(cpm.space, property)',
            c = c,
            grid=false,
            axis=nothing,
            legend=legend,
            framestyle=:box,
            aspect_ratio=:equal,
            size = figureSize,
            xlims=(0.5, rows+0.5),
            ylims=(0.5, columns+0.5),
            clim=(0,colorMax),
            kwargs...
            )

        cellborders!(plotObject, cpm.space)

        cellMovement!(plotObject,cpm, colorMax)

        ModelStep!(cpm)

        plotObject
    end every frameskip

    return gif(anim, file, fps = framerate)

end


function meshscatter(space::CellSpace; legend=false)

    plotObject = plot(
        lims=(0,maximum(size(space))),
        framestyle=:box,
        size = (600, 600)
    )

    xp = [0, 0, 0, 0, 1, 1, 1, 1] .- 0.5;
    yp = [0, 1, 0, 1, 0, 0, 1, 1] .- 0.5;
    zp = [0, 0, 1, 1, 1, 0, 0, 1] .- 0.5;
    connections = [(1,2,4,3),(1,3,5,6),(5,6,7,8),(2,4,8,7),(1,2,7,6),(3,4,8,5)];


    for I in CartesianIndices(space.nodeIDs)
        if !iszero(space.nodeIDs[I])
            (x,y,z) = Tuple(I)
            mesh3d!(plotObject, xp .+ x, yp .+ y, zp .+ z; connections, proj_type = :persp, color=:dodgerblue, linecolor=:black, alpha=1.0, legend=legend)
        end
    end

    return plotObject
end


function recordCPM(
    file::String,
    cpm::CellPotts{3, T, V, U},
    timestamps = 0:3000,
    c = cgrad(:tol_light, rev=true);
    legend=:none,
    framerate=30,
    kwargs...) where {T,V,U}


    anim = @animate for t in timestamps

        plotObject = meshscatter(cpm.space; legend=legend)

        ModelStep!(cpm)

        plotObject
    end every 10

    return gif(anim, file, fps = framerate)

end
