using Revise
using CellularPotts
using Random

Random.seed!(314159)

space = CellSpace(200,25,isPeriodic=false)

initialCellState = CellTable(
    [:Epithelial],
    [300],
    [1]);

positions = [size(space) .÷ 2]

initialCellState = addcellproperty(initialCellState, :positions, positions)


penalties = [
    AdhesionPenalty([0 30;
                    30  0]),
    VolumePenalty([5]),
    MigrationPenalty(100, [100], size(space))
    ]

cpm = CellPotts(space, initialCellState, penalties)



recordCPM("TightSpaces.gif", cpm, timestamps = 0:1000, figureSize = (2000,200))