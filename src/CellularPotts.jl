module CellularPotts

using OffsetArrays, #Allow some arrays to be zero indexed to include medium
      GLMakie, #For plotting and making the GUI (also exports AbstractPlotting)
      Colors, #More color options for the cells (e.g. :Purples)
      StatsBase, Random, #Currently just used for countmap, inverse_rle; shuffle
      SparseArrays, #For creating large adjacency matrices
      LinearAlgebra, #Additiona functionality for arrays
      LightGraphs, #Needed for creating graphs that look like graphDimension
      Metis #Lighting fast method for partitioning graphs


include("BaseFunctions.jl")
include("DivideCells.jl")
include("InitializeCells.jl")
include("MarkovStep.jl")
include("Gui.jl")

export 

#BaseFunctions.jl
      Hamiltonian,
      AdhesionPenalty,
      VolumePenalty,
      ModelParameters,
      NamedGraph,
      CellAttributes,
      CellPotts,
      UpdateConnections!,

#MarkovStep
      MHStepInfo,
      MHStep!,
      MHStep_naive!,

#DivideCells.jl
      CellDivision!,

#Gui.jl
      CellGUI

end
