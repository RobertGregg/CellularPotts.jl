using CellularPotts
using Test
using Documenter
using Graphs
using Random
using BenchmarkTools

#Setting global random seed for tests
Random.seed!(314159)

####################################################
# Adding/Removing cells
####################################################


cpm = CellPotts(
        CellSpace(10,10; periodic=false),
        CellState(names=[:A,:B,:C], volumes=[8,10,15], counts=[1,1,1]),
        [AdhesionPenalty(fill(20,3,3)), VolumePenalty([5,5])]
    )



####################################################
# CellSpace
####################################################

@testset "Spaces or graphs, who's to say?" begin
    g = CellSpace(3,3; periodic=false, diagonal=false)

    @test nv(g) == 9
    @test ne(g) == 12
    @test eltype(collect(edges(g))) == Graphs.SimpleGraphs.SimpleEdge{Int64}
end

@testset "Periodic grids" begin
    g1 = CellSpace(3,3; periodic=false, diagonal=false)
    g2 = CellSpace(3,3; periodic=true, diagonal=false)
    g3 = CellSpace(3,3; periodic=false, diagonal=true)
    g4 = CellSpace(3,3; periodic=true, diagonal=true)

    @test ne(g1) == 12
    @test ne(g2) == 18
    @test ne(g3) == 20
    @test ne(g4) == 36

    @test repr(g1) == "3×3 nonPeriodic 4-Neighbor CellSpace{Int64,2}"    
    @test repr(g2) == "3×3 Periodic 4-Neighbor CellSpace{Int64,2}"    
    @test repr(g3) == "3×3 nonPeriodic 8-Neighbor CellSpace{Int64,2}"    
    @test repr(g4) == "3×3 Periodic 8-Neighbor CellSpace{Int64,2}"    
end

@testset "Neighborhoods" begin
    g1 = CellSpace(3,3; diagonal=false)
    g2 = CellSpace(3,3; diagonal=true)

    @test ne(g1) == 18
    @test ne(g2) == 36
end

####################################################
# CellTables
####################################################

@testset "Lonely cells" begin

    table1 = CellState(:hello, 1,1)
    table2 = CellState([:hello], [1],[1])
    table3 = CellState(names=:hello, volumes=1, counts=1)

    @test parent(table1) == parent(table2) == parent(table3)
    
    @test repr(table1) == """
    ┌─────┬────────┬─────────┬─────────┬─────────┬────────────────┬────────────┬───────────────────┐
    │ Row │  names │ cellIDs │ typeIDs │ volumes │ desiredVolumes │ perimeters │ desiredPerimeters │
    │     │ Symbol │   Int64 │   Int64 │   Int64 │          Int64 │      Int64 │             Int64 │
    ├─────┼────────┼─────────┼─────────┼─────────┼────────────────┼────────────┼───────────────────┤
    │   0 │ Medium │       0 │       0 │       0 │              0 │          0 │                 0 │
    │   1 │  hello │       1 │       1 │       0 │              1 │          0 │                 8 │
    └─────┴────────┴─────────┴─────────┴─────────┴────────────────┴────────────┴───────────────────┘
    """
end

@testset "Testing positions keyword" begin

    table1 = CellState(:hello, 1, 1; positions=(25,25))
    table2 = CellState([:hello], [1], [1]; positions=[(25,25)])
    table3 = CellState(names=:hello, volumes=1, counts=1, positions=(25,25))

    table4 = CellState([:hello], [10], [3]; positions=[(1,1),(2,2),(3,3)])
    table5 = CellState(:hello, 10, 3; positions=[(1,1),(2,2),(3,3)])
    table6 = CellState(names=:hello, volumes=10, counts=3, positions=[(1,1),(2,2),(3,3)])

    @test parent(table1) == parent(table2) == parent(table3)
    @test parent(table4) == parent(table5) == parent(table6)
end

@testset "Adding Cell Properties" begin

    df = CellState(
    names = [:Epithelial, :TCell, :BCell],
    volumes = [500, 100, 100],
    counts = [1, 2, 3],
    P1 = [0,0,0],
    P2 = [0, missing, missing],
    P3 = collect(1:6)
    )

    for i in 1:3
        Pi = Symbol("P$(i)")
        @test hasproperty(df, Pi)
        @test length(getproperty(df, Pi)) == first(size(df))
    end
end

@testset "Add/Delete Cells" begin
    table = CellState(names=:hello, volumes=1, counts=1)

    push!(table,table[1])

    @test countcells(table) == 2

    deleteat!(table,1)

    @test countcells(table) == 1
end

####################################################
# Penalties
####################################################

@testset "Penalty tests" begin

    @test_throws ErrorException("J needs to be symmetric") AdhesionPenalty([1 2; 3 4])
end


####################################################
# Articulation Point Detector
####################################################

@testset "General Fragment Test" begin
    space = CellSpace(5,5; periodic=true, diagonal=true)

    initialCellState = CellState(:Epithelial, 8, 1)

    penalties = [AdhesionPenalty([0 20; 20 0]), VolumePenalty([5])]

    cpm = CellPotts(space, initialCellState, penalties)

    #Rework cpm to have an articulation point
    cpm.space.nodeIDs .= 0
    pos = [[2,2], [2,3], [3,2], [3,3], [3,4], [4,3], [4,4], [5,5]]
    
    for (i,j) in pos
        cpm.space.nodeIDs[i,j] = 1
        cpm.space.nodeTypes[i,j] = 1
    end

    cpm.step.source.node = 24
    cpm.step.source.id = 0
    cpm.step.source.type = 0
    cpm.step.source.neighbors = neighbors(cpm.space,24)

    cpm.step.target.node = 19
    cpm.step.target.id = 1
    cpm.step.target.type = 1
    cpm.step.target.neighbors = neighbors(cpm.space,19)

    @test isfragmented(cpm) == true
end


@testset "Euler Characteristic" begin
    space = CellSpace(5,5; periodic=true, diagonal=false)

    initialCellState = CellState(:Epithelial, 8, 1)

    penalties = [AdhesionPenalty([0 20; 20 0]), VolumePenalty([5])]

    cpm = CellPotts(space, initialCellState, penalties)

    #Rework cpm to have an articulation point
    cpm.space.nodeIDs .= 0
    pos = [[2,2], [2,3], [3,2], [3,3], [3,4], [4,3], [4,4], [5,4]]
    
    for (i,j) in pos
        cpm.space.nodeIDs[i,j] = 1
        cpm.space.nodeTypes[i,j] = 1
    end

    #Test Articulation Point
    cpm.step.source.node = 24
    cpm.step.source.id = 0
    cpm.step.source.type = 0
    cpm.step.source.neighbors = neighbors(cpm.space,24)

    cpm.step.target.node = 19
    cpm.step.target.id = 1
    cpm.step.target.type = 1
    cpm.step.target.neighbors = neighbors(cpm.space,19)

    @test isfragmented(cpm) == true

    #Test not an articulation point
    cpm.step.source.node = 11
    cpm.step.source.neighbors = neighbors(cpm.space,11)

    cpm.step.target.node = 12
    cpm.step.target.neighbors = neighbors(cpm.space,12)

    @test isfragmented(cpm) == false

end

####################################################
# Show method for CPM
####################################################

@testset "Show me the CPM" begin
    
    cpm = CellPotts(
        CellSpace(5,5),
        CellState(names=[:A,:B], volumes=[8,10], counts=[1,1]),
        [AdhesionPenalty([0 20; 20 0]), VolumePenalty([5])]
    )

    @test repr(cpm) == "Cell Potts Model:\nGrid: 5×5\nCell Counts: [A → 1] [B → 1] [Total → 2]\nModel Penalties: Adhesion Volume\nTemperature: 20.0\nSteps: 0"
end

####################################################
# Example Gallery
####################################################

@testset "Hello World" begin

    file = "../docs/src/ExampleGallery/HelloWorld/HelloWorld.jl"

    #Skip expressions that create visualizations
    skipRecord(expr) = first(expr.args) == :record ? nothing : expr

    include(skipRecord, file) #creates the hello world model

    ModelStep!(cpm)

    @test cpm.step.counter == 1
end

@testset "Going 3D" begin

    file = "../docs/src/ExampleGallery/Going3D/Going3D.jl"

    #Skip expressions that create visualizations
    skipRecord(expr) = first(expr.args) == :record ? nothing : expr

    include(skipRecord, file)

    ModelStep!(cpm)

    @test cpm.step.counter == 1
end

@testset "Lets Get Moving" begin

    file = "../docs/src/ExampleGallery/LetsGetMoving/LetsGetMoving.jl"

    #Skip expressions that create visualizations
    skipRecord(expr) = first(expr.args) == :record ? nothing : expr

    include(skipRecord, file) 

    ModelStep!(cpm)

    @test cpm.step.counter == 1
end

@testset "OnPatrol" begin

    file = "../docs/src/ExampleGallery/OnPatrol/OnPatrol.jl"

    #Skip expressions that create visualizations
    skipRecord(expr) = first(expr.args) == :record ? nothing : expr

    include(skipRecord, file) 

    ModelStep!(cpm)

    @test cpm.step.counter > 1 #ModelStep used in script 
end


@testset "Over Here" begin

    file = "../docs/src/ExampleGallery/OverHere/OverHere.jl"

    #Skip expressions that create visualizations
    skipRecord(expr) = first(expr.args) ∈ [:record, :gif, :anim] ? nothing : expr

    include(skipRecord, file) 

    ModelStep!(cpm)

    @test cpm.step.counter == 1 #ModelStep used in script 
end


@testset "Bringing ODEs To Life" begin

    file = "../docs/src/ExampleGallery/BringingODEsToLife/BringingODEsToLife.jl"

    #Skip expressions that create visualizations
    skipRecord(expr) = first(expr.args) ∈ [:record, :gif, :anim, :plot] ? nothing : expr

    include(skipRecord, file) 

    ModelStep!(cpm)

    println(length(last(sol)))

    @test cpm.step.counter > 1 #ModelStep used in script 
end


@testset "Diffusion Inside Cells" begin

    file = "../docs/src/ExampleGallery/DiffusionInsideCells/DiffusionInsideCells.jl"

    #Skip expressions that create visualizations
    skipRecord(expr) = first(expr.args) ∈ [:record, :gif, :anim, :plot, :plot!] ? nothing : expr

    include(skipRecord, file) 

    ModelStep!(cpm)

    @test cpm.step.counter > 1 #ModelStep used in script 
end


@testset "Diffusion Outside Cells" begin

    file = "../docs/src/ExampleGallery/DiffusionOutsideCells/DiffusionOutsideCells.jl"

    #Skip expressions that create visualizations
    skipRecord(expr) = first(expr.args) ∈ [:record, :gif, :anim, :plot, :plot!] ? nothing : expr

    include(skipRecord, file) 

    ModelStep!(cpm)

    @test cpm.step.counter > 1 #ModelStep used in script 
end


@testset "Tight Spaces" begin

    file = "../docs/src/ExampleGallery/TightSpaces/TightSpaces.jl"

    #Skip expressions that create visualizations
    skipRecord(expr) = first(expr.args) ∈ [:record, :gif, :anim, :plot, :plot!] ? nothing : expr

    include(skipRecord, file) 

    ModelStep!(cpm)

    @test cpm.step.counter == 1 #ModelStep used in script 
end


@testset "Travel Time" begin

    file = "../docs/src/ExampleGallery/TravelTime/TravelTime.jl"

    #Skip expressions that create visualizations
    skipRecord(expr) = first(expr.args) ∈ [:record, :visualize, :gif, :anim, :plot, :plot!, scatter] ? nothing : expr

    include(skipRecord, file) 

    ModelStep!(cpm)

    @test cpm.step.counter > 1 #ModelStep used in script 
end

####################################################
# Benchmarks
####################################################

@testset "Allocations (no diagonal)" begin

    cpm = CellPotts(
        CellSpace(50,50),
        CellState(:Epithelial, 25, 10),
        [
            AdhesionPenalty([0 30; 30 0]),
            VolumePenalty([5]),
            PerimeterPenalty([0,10]),
            MigrationPenalty(50, [50], (50,50))
            ]
        )
    
    @test (@ballocated ModelStep!(cpm)) == 0
end


@testset "Allocations (diagonal)" begin

    cpm = CellPotts(
        CellSpace(50,50; diagonal=true),
        CellState(:Epithelial, 25, 10),
        [
            AdhesionPenalty([0 30; 30 0]),
            VolumePenalty([5]),
            PerimeterPenalty([0,10]),
            MigrationPenalty(50, [50], (50,50))
            ]
        )
    
    @test (@ballocated ModelStep!(cpm)) == 0
end




##################################
# Delete me
##################################
# using Plots


# function square(xy)
#     x,y=xy
#     [(NaN, NaN),
#     (x - 0.5, y - 0.5),
#     (x - 0.5, y + 0.5),
#     (x + 0.5, y + 0.5),
#     (x + 0.5, y - 0.5)]
# end

# xdim,ydim = (100,100)
# z1 = [100exp(-((x-xdim/2)^2+(y-ydim/2)^2)/10000) for x in 1:xdim, y in 1:ydim]
# z2 = zeros(xdim,ydim)

# for x in 1:xdim, y in 1:ydim
#     if (x-10)^2 + (y-15)^2 < 25
#         z2[x,y] = 1
#     end


#     if (x-40)^2 + (y-15)^2 < 40
#         z2[x,y] = 1
#     end


#     if (x-60)^2 + (y-60)^2 < 100
#         z2[x,y] = 2
#     end
# end


# plt = contourf(z1,
#     c=:temperaturemap,
#     levels=50,
#     alpha=0.9,
#     linewidth=0,
#     size=(600,600),
#     axis=nothing,
#     aspect_ratio=:equal,
#     framestyle = :box,
#     xlims = (1, xdim),
#     ylims = (1, ydim))


# cellColors = [(:darkblue, 0.5), (:darkgreen,0.5)]

# for (i,fillcolor) in enumerate(cellColors)
#     squares = Shape(vcat((square(I) for I in Iterators.product(1:xdim,1:ydim) if z2[I...] == i)...))
#     plot!(plt,squares, legend=true, fill=fillcolor, linealpha=0, label="cell type $i")
# end


# plt

