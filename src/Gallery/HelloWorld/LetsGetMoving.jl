using CellularPotts


space = CellSpace(200,200)

initialCellState = newCellState(
    [:Epithelial],
    [500],
    [1]);

positions = [(100,100)]

initialCellState = addCellProperty(initialCellState, :positions, positions)

penalties = [
    AdhesionPenalty([0 20;
                    20 100]),
    VolumePenalty([50]),
    PerimeterPenalty([5]),
    MigrationPenalty(20,1,[:Epithelial], space.gridSize)
    ]


cpm = CellPotts(space, initialCellState, penalties);

positionCells!(cpm)

ModelStep!(cpm)

CellGUI(cpm)