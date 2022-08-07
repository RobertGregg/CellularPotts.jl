using Revise
using CellularPotts


space = CellSpace(200,200)

initialCellState = newCellState(
    [:Epithelial, :TCell],
    [500, 400],
    [75, 1])


penalties = [
    AdhesionPenalty([0 20 20;
                    20 20 100;
                    20 100 200]),
    VolumePenalty([30,30]),
    PerimeterPenalty([0,2]),
    MigrationPenalty(500, [0,60], space.gridSize)
    ]


cpm = CellPotts(space, initialCellState, penalties)

positionCellsRandom!(cpm)

CellGUI(cpm)