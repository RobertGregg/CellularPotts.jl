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


@testset "MyPackage" begin
    ... # other tests & testsets
    doctest(MyPackage)
    ... # other tests & testsets
end