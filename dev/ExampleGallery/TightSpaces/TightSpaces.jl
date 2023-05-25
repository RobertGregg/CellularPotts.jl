# # Tight Spaces

# In this example we use a matrix of zeros and ones to denote valid locations that cells may occupy. This example blocks out a square in the middle of the space where cell cannot move to. 

# More complex geometrics (2D and 3D) could be created by modifying the spaceImage defined below.

using CellularPotts
using Random
using Plots

Random.seed!(314159)

#Make a space that looks like a box frame
spaceImage = ones(Int, 100,100)
spaceImage[20:80,20:80] .= 0

space = CellSpace(spaceImage, isPeriodic=false)

initialCellState = CellTable(:Epithelial, 500, 1);

positions = [90,10]

initialCellState = addcellproperty(initialCellState, :positions, positions)


penalties = [
    AdhesionPenalty([0 30;
                    30  0]),
    VolumePenalty([10]),
    MigrationPenalty(100, [100], size(space))
    ]

cpm = CellPotts(space, initialCellState, penalties)





(rows,columns) = size(cpm.space)

anim = @animate for t in 0:1000

    plotObject = heatmap(
        cpm.space.nodeIDs',
        c = cgrad(:tol_light, rev=true),
        grid=false,
        axis=nothing,
        legend=:none,
        framestyle=:box,
        aspect_ratio=:equal,
        size = (600,600),
        xlims=(0.5, rows+0.5),
        ylims=(0.5, columns+0.5)
        )

    plot!(plotObject,[19.5, 19.5],[19.5, 80.5], color=:black)
    plot!(plotObject,[80.5, 80.5],[19.5, 80.5], color=:black)
    plot!(plotObject,[19.5, 80.5],[19.5, 19.5], color=:black)
    plot!(plotObject,[19.5, 80.5],[80.5, 80.5], color=:black)

    cellborders!(plotObject, cpm.space)

    cellMovement!(plotObject,cpm)

    ModelStep!(cpm)

    plotObject
end

gif(anim, "TightSpaces.gif", fps = 30)