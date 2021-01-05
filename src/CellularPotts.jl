module CellularPotts

using GLMakie, AbstractPlotting, Colors

import StatsBase.counts, Printf.@printf

#Some functions needed for all simulations
include("BaseFunctions.jl")

export CellPotts, MHStep!, Edge2Grid

end
