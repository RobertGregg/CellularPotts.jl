# # On Patrol 

# Here we combine ideas from the [HelloWorld.jl](https://robertgregg.github.io/CellularPotts.jl/dev/ExampleGallery/HelloWorld/HelloWorld/) and [LetsGetMoving.jl](https://robertgregg.github.io/CellularPotts.jl/dev/ExampleGallery/LetsGetMoving/LetsGetMoving/) examples to simulate T-Cells patrolling through a dense layer of epithelial cells. 

# Start by loading in the `CellularPotts.jl` package and creating a space where cells can exist.
using CellularPotts
space = CellSpace(200,200)

# Initialize a new `CellTable` with 75 epithelial cells and 5 T-Cells

initialCellState = CellTable(
    [:Epithelial, :TCell],
    [500, 400],
    [75, 5])

# Note that for the `MigrationPenalty` we set the epithelial cell's scaling factor to zero. This effectively removes this penalty from the cell type.

penalties = [
    AdhesionPenalty([0 30 30;
                    30 30 30
                    30 30 30]),
    VolumePenalty([10,10]),
    MigrationPenalty(50, [0,50], size(space))
    ]

# Create a new `CellPotts` model.

cpm = CellPotts(space, initialCellState, penalties)


# Record the simulation
recordCPM("OnPatrol.gif", cpm)