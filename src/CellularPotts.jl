module CellularPotts

using OffsetArrays, #Allow some arrays to be zero indexed to include medium
      GLMakie, #For plotting and making the GUI (also exports AbstractPlotting)
      Colors, #More color options for the cells (e.g. :Purples)
      StatsBase, Random, #Currently just used for countmap, inverse_rle; shuffle
      SparseArrays, #For creating large adjacency matrices
      LinearAlgebra, #Additional functionality for arrays
      Graphs, #Needed for creating graphs that look like graphDimension
      Metis, #lightning fast method for partitioning graphs (only thing in this package that is not Julia)
      TimerOutputs #looking for allocations

import Graphs:
            AbstractSimpleGraph,
            SimpleEdge,
            is_directed,
            edgetype,
            ne,
            nv,
            vertices,
            edges,
            has_vertex,
            has_edge,
            add_edge!

import Base: eltype


#include("BaseFunctions.jl")
include("newBase.jl")
include("CellActions.jl")
#include("InitializeCells.jl")
include("MarkovStep.jl")
include("Gui.jl")

export 

#newBase.jl
      InitialCellState,
      CellSummary,
      Penalty,
      AdhesionPenalty,
      VolumePenalty,
      PerimeterPenalty,
      Parameters,
      initializeCells!,
      CellPotts,
      genAdj,
#CellActions.jl
      CellDivision!,
      CellDeath!,
#MarkovStep.jl
      MHStep!,
#Gui.jl
      CellGUI
end
