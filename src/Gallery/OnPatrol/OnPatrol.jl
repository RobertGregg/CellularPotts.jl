using Revise
using CellularPotts


space = CellSpace(100,100)

initialCellState = newCellState(
    [:Epithelial, :TCell],
    [130, 100],
    [60, 10])


penalties = [
    AdhesionPenalty([0 50 50;
                    50 30 50;
                    50 50 30]),
    VolumePenalty([5,5]),
    PerimeterPenalty([5,5]),
    MigrationPenalty(50, [0,100], space.gridSize)
    ]


cpm = CellPotts(space, initialCellState, penalties)

positionCellsRandom!(cpm)

CellGUI(cpm)