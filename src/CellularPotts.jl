module CellularPotts

using OffsetArrays, #Allow some arrays to be zero indexed to include medium
      GLMakie, #For plotting and making the GUI (also exports AbstractPlotting)
      Colors, #More color options for the cells (e.g. :Purples)
      StatsBase, Random, #Currently just used for countmap, inverse_rle; shuffle
      PrettyTables, #Make nice tables for the cell state
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

import Base: eltype, getproperty, merge

import Tables:
            istable,
            schema,
            Schema,
            columnaccess,
            columns,
            getcolumn,
            columnnames,
            rowaccess,
            rows,
            AbstractColumns,
            AbstractRow


include("Spaces.jl")
include("Penalties.jl")
include("Base.jl")
include("InitializeCells.jl")
include("MarkovStep.jl")
# include("CellActions.jl")
include("Gui_new.jl")

export 

#Spaces
      vonNeumannNeighbors,
      mooreNeighbors,
      CellSpace,
#Base.jl
      newCellState,
      addCellProperty,
      getData,
      CellPotts,
      countCells,
      countCellTypes,
      getTypeID,
#Penalties.jl
      Penalty,
      AdhesionPenalty,
      VolumePenalty,
#InitializeCells.jl
      positionCellsRandom!,
      positionCells!,
#MarkovStep.jl
      MHStep!,
# #CellActions.jl
#       CellDivision!,
#       CellDeath!,
#Gui.jl
      CellGUI
end
