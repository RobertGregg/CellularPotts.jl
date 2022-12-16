```@meta
EditURL = "<unknown>/docs/src/ExampleGallery/OnPatrol/OnPatrol.jl"
```

# On Patrol

Here we combine ideas from the [HelloWorld.jl](https://robertgregg.github.io/CellularPotts.jl/dev/ExampleGallery/HelloWorld/HelloWorld/) and [LetsGetMoving.jl](https://robertgregg.github.io/CellularPotts.jl/dev/ExampleGallery/LetsGetMoving/LetsGetMoving/) examples to simulate T-Cells patrolling through a dense layer of epithelial cells.

Start by loading in the `CellularPotts.jl` package and creating a space where cells can exist.

````julia
using CellularPotts
space = CellSpace(200,200)
````

````
200×200 Periodic 8-Neighbor CellSpace{2,Int64}
````

Initialize a new `CellTable` with 75 epithelial cells and 5 T-Cells

````julia
initialCellState = CellTable(
    [:Epithelial, :TCell],
    [500, 400],
    [75, 5]);
````

Note that for the `MigrationPenalty` we set the epithelial cell's scaling factor to zero. This effectively removes this penalty from the cell type.

````julia
penalties = [
    AdhesionPenalty([0 30 30;
                    30 30 30
                    30 30 30]),
    VolumePenalty([10,10]),
    MigrationPenalty(50, [0,50], size(space))
    ]
````

````
3-element Vector{Penalty}:
 AdhesionPenalty([0 30 30; 30 30 30; 30 30 30])
 VolumePenalty([0, 10, 10])
 MigrationPenalty(50, [0, 0, 50], sparse(Int64[], Int64[], Int64[], 200, 200))
````

Create a new `CellPotts` model.

````julia
cpm = CellPotts(space, initialCellState, penalties)
````

````
Cell Potts Model:
Grid: 200×200
Cell Counts: [Epithelial → 75] [TCell → 5] [Total → 80]
Model Penalties: Adhesion Migration Volume
Temperature: 20.0
Steps: 0
````

Record the simulation

````julia
recordCPM("OnPatrol.gif", cpm)
````

```@raw html
<p style="text-align:center;">
    <img title="" src="https://github.com/RobertGregg/CellularPotts.jl/blob/master/docs/src/ExampleGallery/OnPatrol/OnPatrol.gif?raw=true" width="445">
</p>
```


---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

