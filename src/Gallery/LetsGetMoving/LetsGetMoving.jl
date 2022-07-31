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
    AdhesionPenalty([0 30;
                    30 100]),
    VolumePenalty([8]),
    PerimeterPenalty([8]),
    MigrationPenalty(50, [100], space.gridSize)
    ]


cpm = CellPotts(space, initialCellState, penalties);

cpm.temperature = 10.0

positionCells!(cpm)

CellGUI(cpm)
recordCPM("LetsGetMoving1.gif", cpm)