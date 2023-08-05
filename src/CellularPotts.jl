module CellularPotts

using OffsetArrays, #Allow some arrays to be zero indexed to include medium
      Plots, #Visualization 
      Tables, #Structure for holding cells and their properties
      Colors, #More color options for the cells (e.g. :Purples)
      ColorSchemes, #For custom cell colors
      StatsBase, Random, #Currently just used for countmap, inverse_rle; shuffle
      PrettyTables, #Make nice tables for the cell state
      LinearAlgebra, #Additional functionality for arrays
      SparseArrays, #Improve speed for some large arrays
      Graphs, #Needed for creating graph spaces
      Metis, #lightning fast method for partitioning graphs (only thing in this package that is not Julia)
      Literate #Auto genereate markdown files from julia scripts

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


####################################################
# Global Helper Functions
####################################################
#Kronecker delta function
δ(x, y) = isequal(x,y) ? 1 : 0

#Given a desired cell volume, calculate minimum perimeter on a square lattice
#There are 3 pages of notes behind this equation
#It's related to the minimum perimeter for a polyomino which is 2ceil(2√V)
#TODO Honestly why isn't perimeter the literal perimeter?
estPerimeter(V::Int) = iszero(V) ? 0 : 4ceil(Int,2sqrt(V)-3) + 2ceil(Int,2sqrt(V+1)-4) + 14

#Returns a zero indexed array
offset(x) = OffsetArray(x, Origin(0))

#see https://github.com/JuliaArrays/OffsetArrays.jl/pull/137
#This is enough for this package (b/c we only use 0-indexed vectors)
deleteat!(v::OffsetVector{T, Vector{T}}, i::T) where T<:Integer = deleteat!(parent(v),i-first(v.offsets))

#loop through all unique pairs in v
allpairs(v) = Iterators.filter(i -> isless(i...), Iterators.product(v,v))

####################################################
# Files to include
####################################################

for file in readdir("src/Structures")
      include(joinpath("Structures",file))
end

include("Core.jl")
include("CalculatePenalty.jl")
include("InitializeCells.jl")
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
      addcellproperty,
      addnewcell,
      removecell,
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
      recordCPM,
      #TODO need consistant name capitalization 
      cellborders!,
      cellMovement!,
      plotcpm
end
