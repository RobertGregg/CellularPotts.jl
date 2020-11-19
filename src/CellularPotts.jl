module CellularPotts

using Makie, AbstractPlotting, Colors

import StatsBase.counts, Printf.@printf

#Some functions needed for all simulations
include("BaseFunctions.jl")

export CellPotts, MHStep!

end
