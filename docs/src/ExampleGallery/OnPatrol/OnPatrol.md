```@meta
EditURL = "<unknown>/docs/src/ExampleGallery/OnPatrol/OnPatrol.jl"
```

#On Patrol

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
initialCellState = newCellState(
    [:Epithelial, :TCell],
    [500, 400],
    [75, 5])
````

````
┌────────────┬─────────┬─────────┬─────────┬────────────────┬────────────┬───────────────────┐
│      names │ cellIDs │ typeIDs │ volumes │ desiredVolumes │ perimeters │ desiredPerimeters │
│     Symbol │   Int64 │   Int64 │   Int64 │          Int64 │      Int64 │             Int64 │
├────────────┼─────────┼─────────┼─────────┼────────────────┼────────────┼───────────────────┤
│     Medium │       0 │       0 │       0 │              0 │          0 │                 0 │
│ Epithelial │       1 │       1 │       0 │            500 │          0 │               264 │
│ Epithelial │       2 │       1 │       0 │            500 │          0 │               264 │
│ Epithelial │       3 │       1 │       0 │            500 │          0 │               264 │
│ Epithelial │       4 │       1 │       0 │            500 │          0 │               264 │
│ Epithelial │       5 │       1 │       0 │            500 │          0 │               264 │
│ Epithelial │       6 │       1 │       0 │            500 │          0 │               264 │
│ Epithelial │       7 │       1 │       0 │            500 │          0 │               264 │
│ Epithelial │       8 │       1 │       0 │            500 │          0 │               264 │
│ Epithelial │       9 │       1 │       0 │            500 │          0 │               264 │
│ Epithelial │      10 │       1 │       0 │            500 │          0 │               264 │
│ Epithelial │      11 │       1 │       0 │            500 │          0 │               264 │
│ Epithelial │      12 │       1 │       0 │            500 │          0 │               264 │
│ Epithelial │      13 │       1 │       0 │            500 │          0 │               264 │
│ Epithelial │      14 │       1 │       0 │            500 │          0 │               264 │
│ Epithelial │      15 │       1 │       0 │            500 │          0 │               264 │
│ Epithelial │      16 │       1 │       0 │            500 │          0 │               264 │
│ Epithelial │      17 │       1 │       0 │            500 │          0 │               264 │
│ Epithelial │      18 │       1 │       0 │            500 │          0 │               264 │
│ Epithelial │      19 │       1 │       0 │            500 │          0 │               264 │
│ Epithelial │      20 │       1 │       0 │            500 │          0 │               264 │
│ Epithelial │      21 │       1 │       0 │            500 │          0 │               264 │
│ Epithelial │      22 │       1 │       0 │            500 │          0 │               264 │
│ Epithelial │      23 │       1 │       0 │            500 │          0 │               264 │
│ Epithelial │      24 │       1 │       0 │            500 │          0 │               264 │
│ Epithelial │      25 │       1 │       0 │            500 │          0 │               264 │
│ Epithelial │      26 │       1 │       0 │            500 │          0 │               264 │
│ Epithelial │      27 │       1 │       0 │            500 │          0 │               264 │
│ Epithelial │      28 │       1 │       0 │            500 │          0 │               264 │
│ Epithelial │      29 │       1 │       0 │            500 │          0 │               264 │
│ Epithelial │      30 │       1 │       0 │            500 │          0 │               264 │
│ Epithelial │      31 │       1 │       0 │            500 │          0 │               264 │
│ Epithelial │      32 │       1 │       0 │            500 │          0 │               264 │
│ Epithelial │      33 │       1 │       0 │            500 │          0 │               264 │
│ Epithelial │      34 │       1 │       0 │            500 │          0 │               264 │
│ Epithelial │      35 │       1 │       0 │            500 │          0 │               264 │
│ Epithelial │      36 │       1 │       0 │            500 │          0 │               264 │
│ Epithelial │      37 │       1 │       0 │            500 │          0 │               264 │
│ Epithelial │      38 │       1 │       0 │            500 │          0 │               264 │
│ Epithelial │      39 │       1 │       0 │            500 │          0 │               264 │
│ Epithelial │      40 │       1 │       0 │            500 │          0 │               264 │
│ Epithelial │      41 │       1 │       0 │            500 │          0 │               264 │
│ Epithelial │      42 │       1 │       0 │            500 │          0 │               264 │
│ Epithelial │      43 │       1 │       0 │            500 │          0 │               264 │
│ Epithelial │      44 │       1 │       0 │            500 │          0 │               264 │
│ Epithelial │      45 │       1 │       0 │            500 │          0 │               264 │
│ Epithelial │      46 │       1 │       0 │            500 │          0 │               264 │
│ Epithelial │      47 │       1 │       0 │            500 │          0 │               264 │
│ Epithelial │      48 │       1 │       0 │            500 │          0 │               264 │
│ Epithelial │      49 │       1 │       0 │            500 │          0 │               264 │
│ Epithelial │      50 │       1 │       0 │            500 │          0 │               264 │
│ Epithelial │      51 │       1 │       0 │            500 │          0 │               264 │
│ Epithelial │      52 │       1 │       0 │            500 │          0 │               264 │
│ Epithelial │      53 │       1 │       0 │            500 │          0 │               264 │
│ Epithelial │      54 │       1 │       0 │            500 │          0 │               264 │
│ Epithelial │      55 │       1 │       0 │            500 │          0 │               264 │
│ Epithelial │      56 │       1 │       0 │            500 │          0 │               264 │
│ Epithelial │      57 │       1 │       0 │            500 │          0 │               264 │
│ Epithelial │      58 │       1 │       0 │            500 │          0 │               264 │
│ Epithelial │      59 │       1 │       0 │            500 │          0 │               264 │
│ Epithelial │      60 │       1 │       0 │            500 │          0 │               264 │
│ Epithelial │      61 │       1 │       0 │            500 │          0 │               264 │
│ Epithelial │      62 │       1 │       0 │            500 │          0 │               264 │
│ Epithelial │      63 │       1 │       0 │            500 │          0 │               264 │
│ Epithelial │      64 │       1 │       0 │            500 │          0 │               264 │
│ Epithelial │      65 │       1 │       0 │            500 │          0 │               264 │
│ Epithelial │      66 │       1 │       0 │            500 │          0 │               264 │
│ Epithelial │      67 │       1 │       0 │            500 │          0 │               264 │
│ Epithelial │      68 │       1 │       0 │            500 │          0 │               264 │
│ Epithelial │      69 │       1 │       0 │            500 │          0 │               264 │
│ Epithelial │      70 │       1 │       0 │            500 │          0 │               264 │
│ Epithelial │      71 │       1 │       0 │            500 │          0 │               264 │
│ Epithelial │      72 │       1 │       0 │            500 │          0 │               264 │
│ Epithelial │      73 │       1 │       0 │            500 │          0 │               264 │
│ Epithelial │      74 │       1 │       0 │            500 │          0 │               264 │
│ Epithelial │      75 │       1 │       0 │            500 │          0 │               264 │
│      TCell │      76 │       2 │       0 │            400 │          0 │               236 │
│      TCell │      77 │       2 │       0 │            400 │          0 │               236 │
│      TCell │      78 │       2 │       0 │            400 │          0 │               236 │
│      TCell │      79 │       2 │       0 │            400 │          0 │               236 │
│      TCell │      80 │       2 │       0 │            400 │          0 │               236 │
└────────────┴─────────┴─────────┴─────────┴────────────────┴────────────┴───────────────────┘

````

Note that for the `MigrationPenalty` we set the epithelial cell's scaling factor to zero. This effectively removes this penalty from the cell type.

````julia
penalties = [
    AdhesionPenalty([0 20 20;
                    20 20 30;
                    20 30 50]),
    VolumePenalty([10,10]),
    PerimeterPenalty([0,5]),
    MigrationPenalty(50, [0,100], size(space))
    ]
````

````
4-element Vector{Penalty}:
 AdhesionPenalty([0 20 20; 20 20 30; 20 30 50])
 VolumePenalty([0, 10, 10])
 PerimeterPenalty([0, 0, 5], 0, 0)
 MigrationPenalty(50, [0, 0, 100], sparse(Int64[], Int64[], Int64[], 200, 200))
````

Create a new `CellPotts` model.

````julia
cpm = CellPotts(space, initialCellState, penalties)
````

````
Cell Potts Model:
Grid: 200×200
Cell Counts: [Epithelial → 75] [TCell → 5] [Total → 80]
Model Penalties: Adhesion Migration Perimeter Volume
Temperature: 20.0
Steps: 0
````

Here we did not specify the positions of the cells, so they be randomly added to the space.

````julia
positionCellsRandom!(cpm)
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

