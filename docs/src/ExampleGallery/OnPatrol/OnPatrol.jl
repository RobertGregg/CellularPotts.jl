# # On Patrol 

# Here we combine ideas from the [HelloWorld.jl](https://robertgregg.github.io/CellularPotts.jl/dev/ExampleGallery/HelloWorld/HelloWorld/) and [LetsGetMoving.jl](https://robertgregg.github.io/CellularPotts.jl/dev/ExampleGallery/LetsGetMoving/LetsGetMoving/) examples to simulate T-Cells patrolling through a dense layer of epithelial cells. 

# Start by loading in the `CellularPotts.jl` package and creating a space where cells can exist.
using CellularPotts
space = CellSpace(200,200; diagonal=true)

# Initialize a new `CellState` with 75 epithelial cells and 5 T-Cells

initialCellState = CellState(
    [:Epithelial, :TCell],
    [250, 200],
    [160, 10]);

# Note that for the `MigrationPenalty` we set the epithelial cell's scaling factor to zero. This effectively removes this penalty from the cell type.

penalties = [
    AdhesionPenalty([30 30 30;
                    30 20 30
                    30 30 40]),
    VolumePenalty([30, 30]),
    PerimeterPenalty([0, 5]),
    MigrationPenalty(75, [0, 100], size(space))
    ]

# Create a new `CellPotts` model.

cpm = CellPotts(space, initialCellState, penalties)

# Finally, lets run the model for a few steps to let the initial cell positions to equalibrate.

for i=1:50
    ModelStep!(cpm)
end

# Record the simulation
recordCPM("OnPatrol.gif", cpm;
    property = :nodeTypes, frameskip=10, c=:RdBu_3)