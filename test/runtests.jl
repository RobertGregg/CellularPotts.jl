using Revise
using CellularPotts
using Test
using BenchmarkTools
using GLMakie


parameters = Parameters()

cpm = CellPotts(parameters)

initializeCells!(cpm)

MHStep!(cpm)

#@benchmark MHStep!(cpm)

# @profview fp(cpm)


function fp(cpm)
    for i=1:10_000_000
        MHStep!(cpm)
    end
end



parameters = Parameters(
    gridSize = (10,10),
    state = InitialCellState(name=[:Epidermal], count=[4], volume=[10]),
      penalties = [AdhesionPenalty([0 20; 20 100]), VolumePenalty([5]), PerimeterPenalty([5])]
)

cpm = CellPotts(parameters)

initializeCells!(cpm)


CellGUI(cpm)


