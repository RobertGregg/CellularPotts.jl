using CellularPotts
using BenchmarkTools

space1 = CellSpace(50,50; isPeriodic=true, neighborhood=:moore)
space2 = CellSpace(50,50; isPeriodic=true, neighborhood=:vonNeumann)

initialCellState1 = CellState(:Epithelial, 500, 1)
initialCellState1 = addcellproperty(initialCellState1, :positions,  [(25,25)])

initialCellState2 = CellState(:Epithelial, 500, 1)
initialCellState2 = addcellproperty(initialCellState2, :positions,  [(25,25)])


penalties = [
    AdhesionPenalty([0 20;
                     20 0]),
    VolumePenalty([5])
    ]


cpm1 = CellPotts(space1, initialCellState1, penalties)
cpm2 = CellPotts(space2, initialCellState2, penalties)

cpm2.temperature = 100.0

@btime ModelStep!($cpm1)


@btime ModelStep!($cpm2)