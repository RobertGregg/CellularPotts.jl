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
    VolumePenalty([5])
    ]


cpm = CellPotts(space, initialCellState, penalties);

positionCells!(cpm)

ModelStep!(cpm)

CellGUI(cpm)