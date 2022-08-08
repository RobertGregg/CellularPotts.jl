# #Let's Get Moving

# Many cells have the ability to move within their environment through the contraction of actin filaments. This mechanism leads to cells performing an **"intermittent random walk"** which is characterized by periods of persistent movement, followed by periods of being stationary. 

# Let's see how we can add this kind of behavior to our model,

# Start by loading in the `CellularPotts.jl` package and creating a space where cells can exist.
using CellularPotts


space = CellSpace(200,200)

initialCellState = newCellState(
    [:Epithelial],
    [500],
    [1]);

#positions = [(100,100)]

#initialCellState = addcellproperty(initialCellState, :positions, positions)

penalties = [
    AdhesionPenalty([0 30;
                    30 30]),
    VolumePenalty([10]),
    PerimeterPenalty([10]),
    MigrationPenalty(50, [100], space.gridSize)
    ]


cpm = CellPotts(space, initialCellState, penalties);

cpm.temperature = 10.0

positionCellsRandom!(cpm)

CellGUI(cpm)
recordCPM("LetsGetMoving1.gif", cpm)