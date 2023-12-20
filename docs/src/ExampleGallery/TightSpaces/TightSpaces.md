```@meta
EditURL = "TightSpaces.jl"
```

# Tight Spaces

In this example we use a matrix of zeros and ones to denote valid locations that cells may occupy. This example blocks out a square in the middle of the space where cell cannot move to.

More complex geometrics (2D and 3D) could be created by modifying the spaceImage defined below.

````julia
using CellularPotts
using Plots


#Make a space that looks like a box frame
spaceImage = ones(Int, 100,100)
spaceImage[10:90,10:90] .= 0

space = CellSpace(spaceImage; periodic=false, diagonal=true)

initialCellState = CellState(:Epithelial, 500, 1; positions = (5,5));

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

    plot!(plotObject,[9.5, 9.5],[9.5, 90.5], color=:black)
    plot!(plotObject,[90.5, 90.5],[9.5, 90.5], color=:black)
    plot!(plotObject,[9.5, 90.5],[9.5, 9.5], color=:black)
    plot!(plotObject,[9.5, 90.5],[90.5, 90.5], color=:black)

    cellborders!(plotObject, cpm.space)

    cellmovement!(plotObject,cpm)

    ModelStep!(cpm)

    plotObject
end

gif(anim, "TightSpaces.gif", fps = 30)
````

```@raw html
<p style="text-align:center;">
    <img title="" src="https://github.com/RobertGregg/CellularPotts.jl/blob/master/docs/src/ExampleGallery/TightSpaces/TightSpaces.gif?raw=true" width="445">
</p>
```


---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

