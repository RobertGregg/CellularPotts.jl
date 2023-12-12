using Revise
using CellularPotts
using Test
using Documenter
using Graphs


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

    @test parent(table1) == parent(table2)
end

@testset "Adding Cell Properties" begin

    df = CellState(
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
# Example Gallery
####################################################

@testset "Hello World" begin

    file = "../docs/src/ExampleGallery/HelloWorld/HelloWorld.jl"

    #Skip expressions that create visualizations
    skipRecord(expr) = first(expr.args) == :recordCPM ? nothing : expr

    include(skipRecord, file) #creates the hello world model

    ModelStep!(cpm)

    @test cpm.step.counter == 1
end

@testset "Going 3D" begin

    file = "../docs/src/ExampleGallery/Going3D/Going3D.jl"

    #Skip expressions that create visualizations
    skipRecord(expr) = first(expr.args) == :recordCPM ? nothing : expr

    include(skipRecord, file)

    ModelStep!(cpm)

    @test cpm.step.counter == 1
end

@testset "Lets Get Moving" begin

    file = "../docs/src/ExampleGallery/LetsGetMoving/LetsGetMoving.jl"

    #Skip expressions that create visualizations
    skipRecord(expr) = first(expr.args) == :recordCPM ? nothing : expr

    include(skipRecord, file) 

    ModelStep!(cpm)

    @test cpm.step.counter == 1
end

@testset "OnPatrol" begin

    file = "../docs/src/ExampleGallery/OnPatrol/OnPatrol.jl"

    #Skip expressions that create visualizations
    skipRecord(expr) = first(expr.args) == :recordCPM ? nothing : expr

    include(skipRecord, file) 

    ModelStep!(cpm)

    @test cpm.step.counter > 1 #ModelStep used in script 
end


@testset "Over Here" begin

    file = "../docs/src/ExampleGallery/OverHere/OverHere.jl"

    #Skip expressions that create visualizations
    skipRecord(expr) = first(expr.args) ∈ [:recordCPM, :gif, :anim] ? nothing : expr

    include(skipRecord, file) 

    ModelStep!(cpm)

    @test cpm.step.counter == 1 #ModelStep used in script 
end


@testset "Bringing ODEs To Life" begin

    file = "../docs/src/ExampleGallery/BringingODEsToLife/BringingODEsToLife.jl"

    #Skip expressions that create visualizations
    skipRecord(expr) = first(expr.args) ∈ [:recordCPM, :gif, :anim, :plot] ? nothing : expr

    include(skipRecord, file) 

    ModelStep!(cpm)

    @test cpm.step.counter > 1 #ModelStep used in script 
end


@testset "Diffusion Inside Cells" begin

    file = "../docs/src/ExampleGallery/DiffusionInsideCells/DiffusionInsideCells.jl"

    #Skip expressions that create visualizations
    skipRecord(expr) = first(expr.args) ∈ [:recordCPM, :gif, :anim, :plot, :plot!] ? nothing : expr

    include(skipRecord, file) 

    ModelStep!(cpm)

    @test cpm.step.counter > 1 #ModelStep used in script 
end


@testset "Diffusion Outside Cells" begin

    file = "../docs/src/ExampleGallery/DiffusionOutsideCells/DiffusionOutsideCells.jl"

    #Skip expressions that create visualizations
    skipRecord(expr) = first(expr.args) ∈ [:recordCPM, :gif, :anim, :plot, :plot!] ? nothing : expr

    include(skipRecord, file) 

    ModelStep!(cpm)

    @test cpm.step.counter > 1 #ModelStep used in script 
end


@testset "Tight Spaces" begin

    file = "../docs/src/ExampleGallery/TightSpaces/TightSpaces.jl"

    #Skip expressions that create visualizations
    skipRecord(expr) = first(expr.args) ∈ [:recordCPM, :gif, :anim, :plot, :plot!] ? nothing : expr

    include(skipRecord, file) 

    ModelStep!(cpm)

    @test cpm.step.counter == 1 #ModelStep used in script 
end


@testset "Travel Time" begin

    file = "../docs/src/ExampleGallery/TravelTime/TravelTime.jl"

    #Skip expressions that create visualizations
    skipRecord(expr) = first(expr.args) ∈ [:recordCPM, :gif, :anim, :plot, :plot!, scatter] ? nothing : expr

    include(skipRecord, file) 

    ModelStep!(cpm)

    @test cpm.step.counter > 1 #ModelStep used in script 
end