module CellularPotts

using OffsetArrays, #Allow some arrays to be zero indexed to include medium
      GLMakie, #For plotting and making the GUI (also exports AbstractPlotting)
      Colors, #More color options for the cells (e.g. :Purples)
      Optim, #Optimize the slope of a line dividing a cell
      StatsBase #Currently just used for the `counts()` function

import Printf.@printf
import ShiftedArrays.circshift as circle_shift

#Some functions needed for all simulations
include("BaseFunctions.jl")
include("CellDivision.jl")

export CellPotts, Neighbors, MHStep!, Edge2Grid, Propose!, CellDivide!

end
