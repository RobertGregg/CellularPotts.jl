using CellularPotts
using Test
using Random



@testset "constructor methods" begin
    # We want to create a model in a variety of ways

    #All default values
    CellPotts()

    #Specify one desired volume for all cells
    CellPotts(n = 100, σ = 50, Vd = 20)

    #Random starting grid vs seeded vs rectangular
    CellPotts(initialGrid = :random)
    CellPotts(initialGrid = :seeded)
    CellPotts(initialGrid = :rect)

    #Allowing cell properties
    CellPotts(patrol = true)
    CellPotts(divide = true)

end
