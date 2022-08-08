using Revise
using CellularPotts
using Test
using Documenter
using BenchmarkTools
using Graphs


@testset "Spaces or graphs, who's to say?" begin
    g = CellSpace(3,3;wrapAround=false,cellNeighbors=vonNeumannNeighbors)

    @test nv(g) == 9
    @test ne(g) == 12
    @test eltype(collect(edges(g))) == Graphs.SimpleGraphs.SimpleEdge{Int64}
end


@testset "Adding Properties" begin

    df = newCellState(
    [:Epithelial, :TCell, :BCell],
    [500, 100, 100],
    [1, 1, 1])

    df = addcellproperty(df, :P1, 0)
    df = addcellproperty(df, :P2, [0, 0, 0])
    df = addcellproperty(df, :P3, 0, :Epithelial)
    df = addcellproperty(df, :P4, [0], :Epithelial)
    df = addcellproperty(df, :P5, 0, [:Epithelial])
    df = addcellproperty(df, :P6, 1, [:TCell, :BCell])
    df = addcellproperty(df, :P7, [1,2], [:TCell, :BCell])
    df = addcellproperty(df, :P8, Dict([1,2] .=> [:TCell, :BCell]))


    for i in 1:8
        @test hasproperty(parent(df), Symbol("P$(i)"))
        @test length(getproperty(parent(df), Symbol("P$(i)"))) == length(df)
    end

    @test_throws DimensionMismatch addcellproperty(df, :Pbad, [1,2], [:Epithelial, :TCell, :BCell])
    @test_throws DimensionMismatch addcellproperty(df, :Pbad, [1,2,3], [:TCell, :BCell])
    @test_throws TypeError addcellproperty(df, :P1, [1,2,3], :TCell)
end


# @testset "Docs Testing" begin
#     doctest(CellularPotts)
# end