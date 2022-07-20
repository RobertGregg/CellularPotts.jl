module CellularPotts

using OffsetArrays, #Allow some arrays to be zero indexed to include medium
      GLMakie, #For plotting and making the GUI (also exports AbstractPlotting)
      Colors, #More color options for the cells (e.g. :Purples)
      ColorSchemes, #For custom cell colors
      StatsBase, Random, #Currently just used for countmap, inverse_rle; shuffle
      PrettyTables, #Make nice tables for the cell state
      LinearAlgebra, #Additional functionality for arrays
      SparseArrays, #Improve speed for some large arrays
      Graphs, #Needed for creating graphs that look like graphDimension
      Metis #lightning fast method for partitioning graphs (only thing in this package that is not Julia)

import OffsetArrays: Origin

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

import Base: eltype,
             getproperty,
             merge,
             parent,
             show,
             iterate,
             getindex,
             keys,
             values,
             pairs

include("Spaces.jl")
include("Core.jl")
include("Penalties.jl")
include("InitializeCells.jl")
include("MarkovStep.jl")
include("CellActions.jl")
include("Gui_new.jl")

export 

#Spaces
      vonNeumannNeighbors,
      mooreNeighbors,
      CellSpace,
#Core.jl
      newCellState,
      addCellProperty,
      addNewCell,
      countcells,
      countcelltypes,
      CellPotts,
      getTypeID,
      perimeterLocal,
      Penalty,
      AdhesionPenalty,
      VolumePenalty,
      PerimeterPenalty,
      MigrationPenalty,
#Penalties.jl
      createPenalties,
      addPenalty!,
#InitializeCells.jl
      positionCellsRandom!,
      positionCells!,
#MarkovStep.jl
      MHStep!,
      applyPenalties,
# #CellActions.jl
       CellDivision!,
       CellDeath!,
#Gui.jl
      CellGUI
end
