using Revise
using CellularPotts
using Test
using Random
using Graphs
using GLMakie

parameters = Parameters()

cpm = CellPotts(parameters)

initializeCells!(cpm)

#MHStep!(cpm)

# @profview fp(cpm)


function fp(cpm)
    for i=1:1000
        MHStep!(cpm)
    end
end



CellGUI(cpm)


