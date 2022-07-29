using Revise
using CellularPotts


space = CellSpace(200,200)

initialCellState = newCellState(
    [:Epithelial],
    [500],
    [1]);

positions = [(100,100)]

initialCellState = addCellProperty(initialCellState, :positions, positions)

penalties = [
    AdhesionPenalty([0 10;
                    10 0]),
    VolumePenalty([5]),
    PerimeterPenalty([2]),
    MigrationPenalty([30],[300], space.gridSize)
    ]


cpm = CellPotts(space, initialCellState, penalties);

cpm.temperature = 10.0

positionCells!(cpm)

ModelStep!(cpm)

CellGUI(cpm)