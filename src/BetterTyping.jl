using LightGraphs
import LightGraphs.neighbors
# using Plots
using OffsetArrays
using Metis
using StatsBase #countmap
using Statistics #median
using Printf # @printf
using LinearAlgebra # issymmetric
using Random #for testing


#This is only for grids with all dimensions being the same length
    Base.mod1(x::CartesianIndex, y::Int)  = CartesianIndex(mod1.(x.I,y))
#I struggle so much with types...this works, but is it "correct"?
    Base.mod1(x::CartesianIndex{N}, y::NTuple{N,Int}) where N  = CartesianIndex(mod1.(x.I,y))

####################################################
# Structure for Cell
####################################################

mutable struct Cell
    id::Int
    volume::Int
    type::String
    individualGraph::SimpleGraph{Int}

    Cell(id::Int, volume::Int, type::String) = new(id, volume, type, SimpleGraph())
end


#Would it make sense to treat medium as a cell?

####################################################
# Structure for Grid/Graph
####################################################
abstract type GraphGrid end

mutable struct SquareGraph <: GraphGrid
    dimension::NTuple #the length, width, height, etc of the grid
    boundaryCondition::Union{String, Symbol} #should the grid wrap around?
    network::SimpleGraph{Int} #a network of nodes and edges equivalent to the grid

    function SquareGraph(dimension::NTuple; boundaryCondition = "Periodic")

        #Create a graph with the correct size
        numberOfNodes = prod(dimension)
        numberOfDims = length(dimension)
        network = SimpleGraph(numberOfNodes)

        
        #Index the nodes based on the dimension
        nodes = LinearIndices(dimension)
        nodeIndices = CartesianIndices(nodes)

        #Create "unit vectors" that point in each dimension (e.g. north, south, east, west)
        unitIndices = CartesianIndices(ntuple(x->-1:1,numberOfDims))[:] 
        filter!(x->sum(abs, x.I)==1, unitIndices) #filter out Moore neighbors

        #Add edges to the graph, looping through each node
        for (node, nodeIndex) in zip(nodes,nodeIndices)
            #What nodes are adjacent to our current node
            nodesToConnect = nodes[[mod1(nodeIndex + direction, dimension) for direction in unitIndices]]

            #Loop through and add the edges
            for neighborNode in nodesToConnect
                add_edge!(network,node,neighborNode)
            end

        end

        return new(dimension,boundaryCondition,network)
    end
end

mutable struct HexGraph <: GraphGrid #In the works
end

####################################################
# Structure for Hamiltonian Penalty
####################################################
abstract type Hamiltonian end


mutable struct AdhesionPenalty <: Hamiltonian
    J::OffsetMatrix{Int, Matrix{Int}}

    function AdhesionPenalty(J::Matrix{Int})
        issymmetric(J) ? nothing : error("J needs to be symmetric")
        
        Joff = OffsetArray(J, CartesianIndex(0, 0):CartesianIndex(size(J).-1))
        return new(Joff)
    end
end


mutable struct VolumePenalty <: Hamiltonian
    desiredVolumes::OffsetVector{Int,Vector{Int64}}
    λᵥ::OffsetVector{Int,Vector{Int64}}

    function VolumePenalty(desiredVolumes::Vector{Int}, λᵥ::Vector{Int})
        desiredVolumesOff = OffsetVector([0; desiredVolumes], 0:length(desiredVolumes))
        λᵥOff = OffsetVector([0; λᵥ], 0:length(λᵥ))
        return new(desiredVolumesOff, λᵥOff)
    end
end


####################################################
# Model Structure
####################################################

mutable struct CellPotts{G<:GraphGrid}
    graph::G
    cells::OffsetVector{Cell,Vector{Cell}} #offset so index 0 is Medium
    cellMembership::Matrix{Int} #List of cell IDs assigned to each node
    cellTypes::Dict{String, Int} #List of unique cell types 
    nodeGroups::OffsetVector{Vector{Int}, Vector{Vector{Int}}} #Each index is a cell, each vector contains nodes in that cell 
    disconnectors::Vector{Int} #List of nodes that will fragment cells if changed 
    penalties::Vector{Hamiltonian}
    energy::Int
    temperature::Float64
    stepCount::Vector{Int}


    function CellPotts(graph::GraphGrid, cells::Vector{Cell}, penalties::Vector{H}; temperature::Float64=3.0) where {H<:Hamiltonian}

        #adding Medium to cells provided"
        cells = [Cell(0,0,"Medium"); cells]

        numberOfCells = length(cells) - 1 #don't include Medium
        numberOfNodes = nodeNumber(graph)

        #Map the zero index to Medium
        cells = OffsetVector(cells, 0:numberOfCells)

        #List of unique cell types
        uniqueTypes =  unique(getfield.(cells, :type))
        cellTypes = Dict( uniqueTypes .=> 0:length(uniqueTypes)-1) #0 index for medium

        cellVolumes = [cell.volume for cell in cells]  #throw error if totalVolume > # of nodes
        sum(cellVolumes) > numberOfNodes ? error("The cells are too large to fit onto the grid") : nothing
        
        #Initialize the cells in the graph (will add more methods in the future)
        approxVolume = round(Int, median(cellVolumes))
        (cellMembership, nodeGroups, individualCellGraphs) = GrowCells(graph, numberOfCells, approxVolume)

        #Articulation points are like bridges connecting parts of the graph, if removed the graph will disconnect
        #The individual cell graphs have their own indices g₁ -> [1,2,3...n₁], g₂ -> [1,2,3...n₂], etc
        #This maps those indicies back to cellMembership and reduces it to a single vector
        disconnectors =  mapreduce( (A,x)-> A[x], vcat, nodeGroups, articulation.(individualCellGraphs))

        nodeGroups = OffsetVector(nodeGroups, 0:numberOfCells)
        #Update cell volumes after random partition
        for (cell, nodes, gr) in zip(cells, nodeGroups, individualCellGraphs)
            cell.volume = length(nodes)
            cell.individualGraph = gr
        end

        #Create an instance of the model with zero energy 
        G = typeof(graph)
        CPM = new{G}(graph, cells, cellMembership, cellTypes, nodeGroups, disconnectors, penalties, 0, temperature, [0,0])

        #Update the energy with the given penalties
        CPM.energy = sum([f(CPM) for f in penalties]) #😮

        return CPM
    end
end


####################################################
# Smaller helper function
####################################################

#Given a node, what are its neighbors
neighbors(g::GraphGrid, v::Int) = LightGraphs.neighbors(g.network::AbstractGraph, v::Integer)

#How many cells are there?
nodeNumber(g::T) where {T<:GraphGrid} = nv(g.network)
nodeNumber(CPM::CellPotts{T}) where {T<:GraphGrid} = nodeNumber(CPM.graph)

#Kronecker delta function
#δ(x::T, y::T) where {T<:Number} = isequal(x,y) ? one(T) : zero(T)
δ(x, y) = isequal(x,y) ? 1 : 0

function τ(CPM::CellPotts, σ::Int)
    currentCellType = CPM.cells[σ].type #get the current cell type
    return CPM.cellTypes[currentCellType] #return an index for J
end 

function transferNode!(CPM::CellPotts, sourceNode::Int, sourceCell::Int, targetCell::Int)
    #nodeGroups
    filter!(!isequal(sourceNode),CPM.nodeGroups[sourceCell]) #remove source node from cell
    push!(CPM.nodeGroups[targetCell],sourceNode) #add it to the target cell

    #individualGraphs
    
    return nothing
end

#Sometimes we want to go back to sane indices
noOff(x) = OffsetArrays.no_offset_view(x)

#This is very ugly and maybe one day I'll make it better
#This function takes in the length of a square grid and returns pairs of all adjacent squares (wraps around)
# ╔═══╤═══╗
# ║ 1 │ 3 ║
# ╠═══╪═══╣
# ║ 2 │ 4 ║
# ╚═══╧═══╝
# Edge2Grid(2) = [[2, 1], [4, 3], [1, 2], [3, 4], [2, 1], [4, 3], [2, 4], [1, 3], [4, 2], [3, 1], [2, 4], [1, 3]]

function Edge2Grid(gridSize::Int)
    gridIndices = 1:gridSize^2

    x1 = reverse(reshape(gridIndices,gridSize,gridSize),dims=1)'[:]
    x2 = circshift(x1,gridSize)

    y1 = reverse(reshape(reverse(gridIndices),gridSize,gridSize),dims=2)[:]
    y2 = circshift(y1,gridSize)

    append!(x1,x1[1:gridSize])
    append!(x2,x2[1:gridSize])
    append!(y1,y1[1:gridSize])
    append!(y2,y2[1:gridSize])

    return [[id1,id2] for (id1,id2) in zip([x1;y1],[x2;y2])]
end

####################################################
# Penalty Functions
####################################################

#----Calculate on entire graphgrid----

function (AP::AdhesionPenalty)(CPM::CellPotts)
    #Initialize energy
    H = 0

    #Loop through all edges in network and calculate adhesion adjacencies
    for edge in edges(CPM.graph.network)
        (σᵢ, σⱼ) = (CPM.cellMembership[edge.src], CPM.cellMembership[edge.dst])

        H += AP.J[τ(CPM, σᵢ), τ(CPM, σⱼ)] * (1-δ(σᵢ, σⱼ))
    end

    return H
end

function (VP::VolumePenalty)(CPM::CellPotts)
    #Extract the current cell volumes
    currentVolumes = [cell.volume for cell in CPM.cells] #Offset Vector

    #Create slots to add up penalties for each cell type
    ΔV = similar(VP.λᵥ)
    ΔV .= 0 #Making an Offset zeros vector

    #Loop through cells, see how far they are from a desired volume, and put into appropriate spot
    for (i, cell) in pairs(CPM.cells) #pairs is enumerate for Offset arrays
        cellTypeIdx = CPM.cellTypes[cell.type]
        ΔV[cellTypeIdx] += VP.λᵥ[cellTypeIdx] * (currentVolumes[i] - VP.desiredVolumes[i])^2
    end

    #Calculate the penalty
    return sum(ΔV)
end


#----Calculate on local graphgrid----

function (AP::AdhesionPenalty)(CPM::CellPotts,
                               sourceNode::Int,
                               sourceCell::Int,
                               sourceNodeNeighbors::Vector{Int},
                               targetCell::Int)

    #Just want to look at the local neighborhood
    subGraphIdx = vcat(sourceNode,sourceNodeNeighbors)

    #Create a source and target membership
    sourceMembership = CPM.cellMembership[subGraphIdx]
    targetMembership = copy(sourceMembership)
    targetMembership[1] = targetCell

    #Loop through all edges for the source
    sourceH = 0
    for edge in edges(CPM.graph.network[subGraphIdx])
        (σᵢ, σⱼ) = (sourceMembership[edge.src], sourceMembership[edge.dst])

        sourceH += AP.J[τ(CPM, σᵢ), τ(CPM, σⱼ)] * (1-δ(σᵢ, σⱼ))
    end

    #Loop through all edges for the source
    targetH = 0
    for edge in edges(CPM.graph.network[subGraphIdx])
        (σᵢ, σⱼ) = (targetMembership[edge.src], targetMembership[edge.dst])

        targetH += AP.J[τ(CPM, σᵢ), τ(CPM, σⱼ)] * (1-δ(σᵢ, σⱼ))
    end

    #Calculate the adhesion difference between target and source
    return targetH - sourceH
end

function (VP::VolumePenalty)(CPM::CellPotts,
                             sourceNode::Int,
                             sourceCell::Int,
                             sourceNodeNeighbors::Vector{Int},
                             targetCell::Int)

    #Create a vector of source and target memberships
    sourceTarget = [sourceCell, targetCell]

    #Extract cell types
    λᵥ = [ CPM.cellTypes[cell.type] for cell in CPM.cells[sourceTarget] ]

    #Extract the desired volumes
    desiredVolumes = VP.desiredVolumes[sourceTarget]

    #Extract the current volumes
    currentVolumes = [ cell.volume for cell in CPM.cells[sourceTarget] ]

    #Calculate the target volumes
    targetVolumes = currentVolumes .+ [-1, 1]
    
    #Return the delta (Target - Source)
    return sum(@. λᵥ*(targetVolumes - desiredVolumes)^2) - sum(@. λᵥ*(currentVolumes - desiredVolumes)^2)
end

####################################################
# Metropolis–Hasting Step
####################################################

function MHStep!(CPM::CellPotts{T}) where {T<:GraphGrid}

    #Pick a random location on the graph
    sourceNode = rand(1:nodeNumber(CPM))
    #What cell does it belong to?
    sourceCell = CPM.cellMembership[sourceNode]

    #Get all of the unique cell IDs neighboring this Node
    sourceNodeNeighbors = neighbors(CPM.graph, sourceNode)
    possibleCellTargets = unique(CPM.cellMembership[sourceNodeNeighbors])

    #Some checks before attempting a flip
    if all(possibleCellTargets .== sourceCell) #In the middle of a cell
        CPM.stepCount[1] += 1
        return nothing

    elseif sourceNode ∈ CPM.disconnectors #will fragment the cell
        CPM.stepCount[1] += 1
        return nothing
    end     

    #Choose a target
    targetCell = rand(possibleCellTargets)

    #Calculate the change in energy when source node is modified
    ΔH =  sum([f(CPM, sourceNode, sourceCell, sourceNodeNeighbors, targetCell) for f in CPM.penalties])

    #Calculate an acceptance ratio
    acceptRatio = min(1.0,exp(-ΔH/CPM.temperature))

    if rand() < acceptRatio #If the acceptance ratio is large enough

        #Update cell volumes
        CPM.cells[sourceCell].volume -= 1
        CPM.cells[targetCell].volume += 1

        #Update cell memberships
        CPM.cellMembership[sourceNode] = targetCell

        #Update energy
        CPM.energy += ΔH

        #Update graphs and node collects
        transferNode!(CPM::CellPotts, sourceNode::Int, sourceCell::Int, targetCell::Int)
        CPM.cells[sourceCell].individualGraph = CPM.graph.network[ CPM.nodeGroups[sourceCell] ]
        CPM.cells[targetCell].individualGraph = CPM.graph.network[ CPM.nodeGroups[targetCell] ]

        #Update disconnectors
        bridgePoints = [articulation(cell.individualGraph) for cell in CPM.cells]
        CPM.disconnectors =  mapreduce( (A,x)-> A[x], vcat, noOff(CPM.nodeGroups), noOff(bridgePoints)) #mapreduce hates offset
        
        #Keep track of steps
        CPM.stepCount[2] += 1
    end
    return nothing
end



function GrowCells(gr::GraphGrid, numberOfCells::Int, approxVolume::Int)

    #Center the cells on the grid
    gridCenter = gr.dimension .÷ 2

    #Create a box in the center of the grid
    totalNodes = numberOfCells * approxVolume #number of grid square in box
    boxDim = ceil(Int,sqrt(totalNodes)) #side length of box

    #Indices for the box
    boxIdx = CartesianIndices((range(gridCenter[1]-boxDim÷2, step=1, length=boxDim),
                            range(gridCenter[2]-boxDim÷2, step=1, length=boxDim) ))

    #Equivalent Linear Indices in the box
    networkIdx = vec(LinearIndices(gr.dimension)[boxIdx])

    #Grid same size as graph
    cellMembership = zeros(Int,gr.dimension)

    #Partition the box by the number of cells
    #The recursive algorithm seems to partition more evenly 
    cellMembership[boxIdx] = Metis.partition(gr.network[networkIdx], numberOfCells, alg=:RECURSIVE)

    
    #Vector of Vectors (1st has all nodes in first cell, etc)
    nodeGroups = [findall(isequal(i),vec(cellMembership)) for i in 0:numberOfCells] #includes Medium

    #Generate subgraphs for each cell
    individualCellGraphs = [gr.network[i] for i in nodeGroups]
    #Check if all cells are connected
    allConnect = sum(is_connected.(individualCellGraphs)) - 1 #remove medium
    allConnect ≠ numberOfCells ? error("some cells are disconnected, try rerunning or use a different cell initialization") : nothing

    return cellMembership, nodeGroups, individualCellGraphs
end






####################################################
# Override Base.show for each struct
####################################################
function Base.show(io::IO, CPM::CellPotts) 
    println("Cell Potts Model:")
    #Grid
    println("$(size(CPM.cellMembership)[1])×$(size(CPM.cellMembership)[2]) $(typeof(CPM.graph))")

    #Cells and types
    print("Cell Counts:")
    allCellTypes = [cell.type for cell in CPM.cells[1:end]]
    for (key, value) in countmap(allCellTypes) #remove medium
        print(" [$(key) → $(value)]")
    end

    if length(unique(allCellTypes)) > 1
        print(" [Total → $(length(CPM.cells))]")
    else
        print("\n")
    end

    print("Model Penalties:")
    for p in typeof.(CPM.penalties)
        print(" $(p)")
    end
    print("\n")
    println("Current Energy: ", CPM.energy)
    println("Grid Temperature: ", CPM.temperature)
end

####################################################
# Testing
####################################################

#Create a graph
gr = SquareGraph((210,210))

#Create a collection of cells
numberOfCells = 10
approxVolume = 2500
cells = [Cell(i,approxVolume,"TCell") for i in 1:numberOfCells]

#Penalties to add to the model
AP = AdhesionPenalty([0 1; 1 1])
VP = VolumePenalty(fill(approxVolume,numberOfCells), [1])

penalties = [AP, VP]

CPM = CellPotts(gr, cells, penalties; temperature = 15.0)

# for i=1:100_000
#     MHStep!(CPM)
#     i % 10_000 == 0 ? println(i) : nothing
# end

#heatmap(CPM.cellMembership)

using BenchmarkTools
Random.seed!(73);
@benchmark MHStep!(CPM)

Random.seed!(73);
@time MHStep!(CPM)
Random.seed!(73);
@time MHStep!(CPM)

# using Profile
# Random.seed!(73);
# @profview MHStep!(CPM)


@code_warntype MHStep!(CPM)

Random.seed!(72);
@time MHStep!(CPM)

a = zeros(100_000)
for i in 1:length(a)
    tmp = @timed MHStep!(CPM)
    a[i] = tmp[:time]
end

histogram(log10.(a))

using Plots




