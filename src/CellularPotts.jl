module CellularPotts

using OffsetArrays, #Allow some arrays to be zero indexed to include medium
      StatsBase, #Sampling weighted vectors, rle_inverse, countmap
      Plots, #Visualization 
      Colors, #More color options for the cells (e.g. :Purples)
      ColorSchemes, #For custom cell colors
      PrettyTables, #Make nice tables for the cell state
      Graphs, #Needed for creating graph spaces
      Metis, #lightning fast method for partitioning graphs (only thing in this package that is not  pure Julia)
      Literate, #Auto genereate markdown files from julia scripts
      LinearAlgebra, #Additional functionality for arrays
      Random, #Currently just used for shuffle
      SparseArrays #Improve speed for some large arrays

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
            add_edge!,
            SimpleGraphs.SimpleEdgeIter,
            SimpleGraphs.fadj,
            SimpleGraphs.badj

import Base: eltype,
             getproperty,
             merge,
             parent,
             show,
             iterate,
             getindex,
             length,
             size,
             deleteat!


####################################################
# Files to include
####################################################

#Global helper functions
include("Miscellaneous.jl")

#Structs
for (root, dirs, files) in walkdir(joinpath(@__DIR__, "Structures"))
      for file in files
            include(joinpath(root, file))
      end
end

include("Core.jl")
include("CalculatePenalty.jl")
include("InitializeCells.jl")
include("Fragmentation.jl")
include("HistoryUpdater.jl")
include("MarkovStepper.jl")
include("CellActions.jl")
include("Visualization.jl")

export 

#Articulation.jl
      Articulation,
      isfragmented,
#CellSpace.jl
      CellSpace,
#CellState.jl
      CellState,
#Penalty.jl
      Penalty,
      AdhesionPenalty,
      VolumePenalty,
      PerimeterPenalty,
      MigrationPenalty,
      ChemoTaxisPenalty,
#Core.jl
      CellPotts,
      countcells,
      countcelltypes,
      cellneighbors,
#Penalties.jl
      addPenalty!,
      perimeterLocal,
#InitializeCells.jl
      positionCellsRandom!,
      positionCells!,
      growcells,
#MarkovStep.jl
      MHStep!,
      calculateΔH,
      ModelStep!,
      updateMHStep!,
      updateModelStep!,
#CellActions.jl
       CellDivision!,
       CellDeath!,
#Visualization.jl
      cellborders,
      visualize,
      visualize!,
      record
end
