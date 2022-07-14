using Revise
using CellularPotts
using Test
using BenchmarkTools
using Graphs


@testset "Spaces or Graphs?" begin
    g = CellSpace(3,3;wrapAround=false,cellNeighbors=vonNeumannNeighbors)

    @test nv(g) == 9
    @test ne(g) == 12
    @test eltype(collect(edges(g))) == Graphs.SimpleGraphs.SimpleEdge{Int64}
end


df = newCellState(
    [:Macrophage, :TCells, :BCells],
    [10, 12, 15],
    [9, 10, 11])


addCellProperty!(df, :TNF, 0.0,[:Macrophage, :TCells])