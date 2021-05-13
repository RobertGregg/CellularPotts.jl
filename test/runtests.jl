using CellularPotts
using Test
using Random

@testset "constructor methods" begin
    # We want to create a model in a variety of ways
end




M = ModelParameters(
    graphDimension = (50,50,50),
    isPeriodic = true,
    cellTypes = ["TCell", "BCell"],
    cellCounts = [100, 50],
    cellVolumes = [100, 80],
    penalties = [AdhesionPenalty([0 1 1; 1 1 1; 1 1 1]),
                 VolumePenalty(fill(900,1500),[1,1])],
    temperature = 20.0)
CPM = CellPotts(M)