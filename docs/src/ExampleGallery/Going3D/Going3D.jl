# # Going 3D
using CellularPotts


space = CellSpace(30,30,30; isPeriodic=false)


initialCellState = CellState(:Epithelial, 1000, 1)


positions = [size(space).÷2]

initialCellState = addcellproperty(initialCellState, :positions, positions)

penalties = [
    AdhesionPenalty([0 30;
                     30 0]),
    VolumePenalty([5])
    ]


cpm = CellPotts(space, initialCellState, penalties)

recordCPM("Going3D.gif", cpm)