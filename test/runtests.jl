using Revise
using CellularPotts
using Test
using Documenter
using BenchmarkTools
using Graphs


####################################################
# CellSpace
####################################################

@testset "Spaces or graphs, who's to say?" begin
    g = CellSpace(3,3; isPeriodic=false, neighborhood=:vonNeumann)

    @test nv(g) == 9
    @test ne(g) == 12
    @test eltype(collect(edges(g))) == Graphs.SimpleGraphs.SimpleEdge{Int64}
end

@testset "Periodic grids" begin
    g1 = CellSpace(3,3; isPeriodic=false)
    g2 = CellSpace(3,3; isPeriodic=true)

    @test ne(g1) == 20
    @test ne(g2) == 36
end

@testset "Neighborhoods" begin
    g1 = CellSpace(3,3; neighborhood=:vonNeumann)
    g2 = CellSpace(3,3; neighborhood=:moore)

    @test ne(g1) == 18
    @test ne(g2) == 36
end

####################################################
# CellTables
####################################################

@testset "Lonely cells" begin

    table1 = CellTable(:hello, 1,1)
    table2 = CellTable([:hello], [1],[1])

    @test parent(table1) == parent(table2)
end

@testset "Adding Cell Properties" begin

    df = CellTable(
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


####################################################
# Penalties
####################################################

@testset "Penalty tests" begin

    @test_throws ErrorException("J needs to be symmetric") AdhesionPenalty([1 2; 3 4])
end