```@meta
EditURL = "OverHere.jl"
```

# Over Here!

In this example we'll demonstrate how to create a static gradient field that attracts cells to its center.

We'll start by importing packages and specifying dimensions for the space cells will occupy. Here those dimensions are declared constants because they will be used repeatedly.

````julia
using CellularPotts, Plots

const xdim = 200
const ydim = 200
space = CellSpace(xdim, ydim, periodic=false, diagonal=true)
````

````
200×200 nonPeriodic 8-Neighbor CellSpace{Int64,2}
````

Initialize 10 cells randomly positioned.

````julia
initialCellState = CellState(:TCells, 200, 10)
````

````
┌────────┬─────────┬─────────┬─────────┬────────────────┬────────────┬───────────────────┐
│  names │ cellIDs │ typeIDs │ volumes │ desiredVolumes │ perimeters │ desiredPerimeters │
│ Symbol │   Int64 │   Int64 │   Int64 │          Int64 │      Int64 │             Int64 │
├────────┼─────────┼─────────┼─────────┼────────────────┼────────────┼───────────────────┤
│ Medium │       0 │       0 │       0 │              0 │          0 │                 0 │
│ TCells │       1 │       1 │       0 │            200 │          0 │               168 │
│ TCells │       2 │       1 │       0 │            200 │          0 │               168 │
│ TCells │       3 │       1 │       0 │            200 │          0 │               168 │
│ TCells │       4 │       1 │       0 │            200 │          0 │               168 │
│ TCells │       5 │       1 │       0 │            200 │          0 │               168 │
│ TCells │       6 │       1 │       0 │            200 │          0 │               168 │
│ TCells │       7 │       1 │       0 │            200 │          0 │               168 │
│ TCells │       8 │       1 │       0 │            200 │          0 │               168 │
│ TCells │       9 │       1 │       0 │            200 │          0 │               168 │
│ TCells │      10 │       1 │       0 │            200 │          0 │               168 │
└────────┴─────────┴─────────┴─────────┴────────────────┴────────────┴───────────────────┘

````

Now to create the static concentration gradient. Higher values will be positioned in the center of the space and expoentially drop off as it moves towards the edges.

````julia
species = [100exp(-((x-xdim/2)^2+(y-ydim/2)^2)/10000) for x in 1:xdim, y in 1:ydim];
````

The new penalty added here is the `ChemoTaxisPenalty`. This penalty takes two arguments, the strength of the chemoattraction and the concentration gradient. Here we create a strong gradient by providing a large scaling factor.

````julia
penalties = [
    AdhesionPenalty([30 30; 30 30]),
    VolumePenalty([10]),
    ChemoTaxisPenalty([50], species)
];
````

Now we can take these three objects and create a Cellular Potts Model object.

````julia
cpm = CellPotts(space, initialCellState, penalties)
````

````
Cell Potts Model:
Grid: 200×200
Cell Counts: [TCells → 10]
Model Penalties: Adhesion Volume ChemoTaxis
Temperature: 20.0
Steps: 0
````

To scale the color gradients properly, we need to know the number of cells in the model

````julia
numCells = countcells(cpm);
````

Normalize the concentration gradient to the number of cells for visualization

````julia
normSpecies = (numCells/maximum(species)) .* species

anim = @animate for t in 1:2000

    plotObject = contourf(normSpecies,
    c=:temperaturemap,
    levels=50,
    alpha=0.2,
    linewidth=0,
    grid=false,
    axis=nothing,
    legend=:none,
    framestyle=:box,
    aspect_ratio=:equal,
    size = (600,600),
    xlims=(0.5, xdim+0.5),
    ylims=(0.5, ydim+0.5))


    heatmap!(plotObject,
    cpm.space.nodeIDs',
    c = cgrad([RGBA(1,1,1,0), RGBA(0,0,0,0.3)], [0.5/numCells], categorical=true),
    clim = (0,numCells)
    )

    cellborders!(plotObject, cpm.space)

    ModelStep!(cpm)

    plotObject

end every 10


gif(anim, "OverHere.gif", fps = 60)
````

```@raw html
<p style="text-align:center;">
    <img title="" src="https://github.com/RobertGregg/CellularPotts.jl/blob/master/docs/src/ExampleGallery/OverHere/OverHere.gif?raw=true" width="445">
</p>
```


---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

