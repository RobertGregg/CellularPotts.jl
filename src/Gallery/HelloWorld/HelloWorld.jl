using Revise
using CellularPotts


space = CellSpace(50,50)

initialCellState = newCellState(
    [:Epithelial],
    [500],
    [1])

positions = [(25,25)]

initialCellState = addcellproperty(initialCellState, :positions, positions)

penalties = [
    AdhesionPenalty([0 20;
                     20 0]),
    VolumePenalty([5])
    ]


cpm = CellPotts(space, initialCellState, penalties)

positionCells!(cpm)

recordCPM("HelloWorld.gif", cpm)