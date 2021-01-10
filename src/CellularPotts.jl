module CellularPotts

using OffsetArrays, GLMakie, AbstractPlotting, Colors
using StatsBase
import Printf.@printf

#Some functions needed for all simulations
include("BaseFunctions.jl")

export CellPotts, MHStep!, Edge2Grid, Propose!

end
