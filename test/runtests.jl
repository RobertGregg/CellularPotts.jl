using CellularPotts
using Test
using Random

@testset "constructor methods" begin
    # We want to create a model in a variety of ways
end




M = ModelParameters()
CPM = CellPotts(M)

MHStep!(CPM)

CellGUI(CPM)