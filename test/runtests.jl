using CellularPotts
using Test
using Random

@testset "constructor methods" begin
    # We want to create a model in a variety of ways
end



# Default 2D
M = ModelParameters()

CPM = CellPotts(M)

MHStep!(CPM)

using BenchmarkTools

@benchmark MHStep!(CPM)


#3D
M = ModelParameters(
    graphDimension = (50,50,50),
    cellCounts = [50],
    cellVolumes = [1000],
    penalties = [AdhesionPenalty([0 50; 50 100]), VolumePenalty(fill(200,50),[5])])

CellGUI(CPM)
