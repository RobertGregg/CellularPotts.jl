using LightGraphs, LinearAlgebra, SparseArrays, Metis
using StatsBase, Random #countmap, inverse_rle; shuffle
using OffsetArrays

#= Notes

Questions
    - what is the difference b/w Rect3D and FRect3D?
    - Issue with visual when dimensions are not equal?
        - An ⊗ I(m) + I(n) ⊗ Am  ≠  Am ⊗ I(n) + I(m) ⊗ An when n≠m

Comments
    - metagraph saves attributes as Dict{Symbol, Any} which leads to a lot of type instability
        - The upside is you can dynamically add any number of attributes
    - Medium could not be connected if not periodic...
        - set a warning for user to try expanding the size of grid
        - basically cells cannot be allowed on the border
    - If you see an offset array, its just to give medium a zero index

Improvements
    - allow additional parameters to nodes (maybe input as a named tuple?)
    - allow cells of the same type to be different sizes
    - could get a big speed improvement if you don't loop through all cells to update articulation points
    - VP having desired volumes makes cell division difficult (type unstable), replace with CPM.M.cellVolumes
    - For 3D gui, don't recreate voxels every iteration, just set color to clear

=#

####################################################
# Helper Functions
####################################################

⊗(A,B) = kron(A,B) #just makes things look nice

#Kronecker delta function
#δ(x::T, y::T) where {T<:Number} = isequal(x,y) ? one(T) : zero(T)
δ(x, y) = isequal(x,y) ? 1 : 0



#circulant and offDiags generate blocks for the adjacency matricies
#See https://stackoverflow.com/a/45958661

function circulant!(A,n)

    vec = spzeros(Int, n)
    vec[[2,n]] .= 1 #put values into the 2nd and last entry

    for i in 1:n
        A[i,:] = vec
        circshift!(vec, vec, 1) #shift the row over by 1
    end

    return A
end

#Fills the diagonals above and below main diagonal
function offDiags!(A)
    A[diagind(A,1)] .= 1
    A[diagind(A,-1)] .= 1
    return A
end

#Generates a square (n,n) sparse adjacency matrix to convert into a network
function genAdj(n::Int, isPeriodic::Bool)
    A = spzeros(Int, n,n)
    
    if isPeriodic #Does the graph wrap around?
        circulant!(A,n)
    else
        offDiags!(A)
    end

    return A ⊗ I(n) + I(n) ⊗ A
end

#Generates a rectangular (m,n) sparse adjacency matrix to convert into a network
function genAdj(m::Int, n::Int, isPeriodic::Bool)
    Am = spzeros(Int, m, m)
    An = spzeros(Int, n, n)

    if isPeriodic #Does the graph wrap around?
        circulant!(Am,m)
        circulant!(An,n)
    else
        offDiags!(Am)
        offDiags!(An)
    end

    return An ⊗ I(m) + I(n) ⊗ Am
end

#Generates a 3D (l,m,n) sparse adjacency matrix to convert into a network
function genAdj(n::Int, m::Int, l::Int, isPeriodic::Bool)
    Al = spzeros(Int, l,l)
    Am = spzeros(Int, m,m)
    An = spzeros(Int, n,n)
    
    if isPeriodic #Does the graph wrap around?
        circulant!(Al,l)
        circulant!(Am,m)
        circulant!(An,n)
    else
        offDiags!(Al)
        offDiags!(Am)
        offDiags!(An)
    end

    return I(l) ⊗ (Am ⊗ I(n) + I(m) ⊗ An) + Al ⊗ I(m*n) #This took a while to figure out...
end

####################################################
# Hamiltonian Penalties
####################################################
abstract type Hamiltonian end


mutable struct AdhesionPenalty <: Hamiltonian
    J::OffsetMatrix{Int, Matrix{Int}} #J[n,m] gives the adhesion penality for cells with types n and m

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
# Structure for user to provide information
####################################################
Base.@kwdef struct ModelParameters{N} #N is graph dimension (makes types stable)
    graphDimension::NTuple{N, Int} = (100,100) #What are the dimensions of the space?
    isPeriodic::Bool = true #Do you want the space to wrap around?
    cellTypes::Vector{String} = ["TCell"] #What kinds of cells do you want?
    cellCounts::Vector{Int} = [10] #How many of each cell do you want?
    cellVolumes::Vector{Int} = [200] #Desired size for each cell? (need same type cells with different sizes?)
    penalties::Vector{Hamiltonian} = [AdhesionPenalty([0 1; 1 6]), VolumePenalty(fill(200,10),[1])] #What penalties are added to the model?
    temperature::Float64 = 3.0 #Temperature to weight site flips (higher increases chance of accepting energy increasing flips)
end


####################################################
# Structure for graph information
####################################################
mutable struct NamedGraph
    network::SimpleGraph{Int} #a network of nodes and edges equivalent to the grid
    #Attributes for the network (length == number of nodes)
    σ::Vector{Int}
    τ::Vector{String}
    isArticulation::BitVector

    #Two inner contructors, the first assume you have all the individual fields
    NamedGraph(network::SimpleGraph{Int}, σ::Vector{Int}, τ::Vector{String}, isArticulation::BitVector) = new(network, σ, τ, isArticulation)

    function NamedGraph(gr::SimpleGraph{Int}, cellMembership::Array{Int})

        #Cell ID for each node
        σ = vec(cellMembership)

        #Initialize a list of cell types proportional to the number of that type
        #Shuffle is to make use all the same times are not next to each other
        cellAssignments = shuffle(inverse_rle(M.cellTypes, M.cellCounts)) #inverse_rle(["a","b"], [2,3]) = ["a","a","b","b","b"] 

        #Each node has a cell type
        τ = [σᵢ == 0 ? "Medium" : cellAssignments[σᵢ] for σᵢ in σ]

        #Additionally, each node can be a articulation point (removing will disconnect the cell)
        isArticulation = falses(length(σ))

        graph = new(gr, σ, τ, isArticulation)

        #Update the articulation points
        UpdateConnections!(graph; checkConnect=true)

        return graph
    end
end

####################################################
# Structure to hold cell level attributes
####################################################
mutable struct CellAttributes
    #Vectors are offset to include medium (medium gets index 0, cell 1 gets index 1, etc.)
    ids::OffsetVector{Int,Vector{Int}} #vector of cell IDs
    volumes::OffsetVector{Int,Vector{Int}} #vector of cell volumes 
    types::OffsetVector{String,Vector{String}} #given a cell index, output it's type
    typeMap::Dict{String, Int} #mapping cell types (e.g. "Medium") to an index (e.g. 0)

    function CellAttributes(graph::NamedGraph, M::ModelParameters)
        #How many cells are there?
        totalCells = sum(M.cellCounts)

        ids = OffsetVector(collect(0:totalCells), 0:totalCells)

        volumes = OffsetVector([sum(x->x==i, graph.σ) for i in 0:totalCells], 0:totalCells)

        types = OffsetVector(["Medium"; inverse_rle(M.cellTypes, M.cellCounts)], 0:totalCells)

        typeMap = Dict( M.cellTypes .=> 1:length(M.cellTypes) )
        typeMap["Medium"] = 0

        return new(ids, volumes, types, typeMap)
    end
end


####################################################
# Model Structure
####################################################

mutable struct CellPotts{N}
    M::ModelParameters{N} #input parameters from the user
    cell::CellAttributes #summary of cell attributes
    graph::NamedGraph #node properties: cell membership, cell type, bridge status
    energy::Int #Total penality energy across graph
    visual::Array{Int,N} #An array of cell memberships for plotting
    stepCounter::Int #counts the number of MHSteps performed

    function CellPotts(M::ModelParameters{N}) where N

        #Generate an adjacency matrix
        adjMat = genAdj(M.graphDimension..., M.isPeriodic)

        #create a graph based on the adjacency matrix
        gr = SimpleGraph(adjMat)

        #Choose a cell initialization method
        cellMembership = GrowCells(gr, M)

        #Add attributes to the nodes
        graph = NamedGraph(gr, cellMembership)

        #create a summary of the cells
        cell = CellAttributes(graph, M)

        #Create an instance of the model with zero energy 
        CPM = new{N}(M, cell, graph, 0, cellMembership, 0)

        #Calculate the energy from the given penalties
        #Update the energy with the given penalties
        CPM.energy = sum([f(CPM) for f in M.penalties]) #😮

        return CPM
    end
end

####################################################
# Variables for Markov Step 
####################################################

Base.@kwdef mutable struct MHStepInfo
    sourceNode::Int = 1
    sourceCell::Int = 1
    sourceNodeNeighbors::Vector{Int} = [1]
    targetCell::Int = 1
end

####################################################
# Penalty Functions
####################################################

#----Calculate on entire graphgrid----

(AP::AdhesionPenalty)(CPM::CellPotts) = AP(CPM.graph, CPM.cell.typeMap)

#This code gets repeated a lot so I made it a separate function
function (AP::AdhesionPenalty)(graph::NamedGraph, typeMap::Dict{String, Int64})

    #Initialize energy
    H = 0

    #Loop through all edges in network and calculate adhesion adjacencies
    for edge in edges(graph.network)
        
        #Given a node index, get the cell id and type
        (σᵢ, σⱼ) = graph.σ[ [edge.src, edge.dst] ]
        (typeᵢ, typeⱼ) = graph.τ[ [edge.src, edge.dst] ]

        #Convert the type (string) to an index for J
        (τᵢ, τⱼ) = ( typeMap[typeᵢ], typeMap[typeⱼ] )

        H += AP.J[τᵢ, τⱼ] * (1-δ(σᵢ, σⱼ))
    end

    return H
end


function (VP::VolumePenalty)(CPM::CellPotts)
    #Create slots to add up penalties for each cell type
    ΔV = similar(VP.λᵥ)
    ΔV .= 0 #Making an Offset zeros vector

    #Loop through cells, see how far they are from a desired volume, and put into appropriate spot
    for id in CPM.cell.ids 
        cellTypeIdx = CPM.cell.typeMap[ CPM.cell.types[id] ] # cell id (21) → cell type ("Medium") → cell type index (0)
        ΔV[cellTypeIdx] += VP.λᵥ[cellTypeIdx] * (CPM.cell.volumes[id] - VP.desiredVolumes[id])^2
    end

    #Calculate the penalty
    return sum(ΔV)
end


#----Calculate on local graphgrid----

function (AP::AdhesionPenalty)(CPM::CellPotts, stepInfo::MHStepInfo)
        
    #Just want to look at the local neighborhood
    subGraphIdx = vcat(stepInfo.sourceNode, stepInfo.sourceNodeNeighbors)

    subgraph = CPM.graph[subGraphIdx] #uses the custom getindex to subset all fields

    sourceH = AP(subgraph, CPM.cell.typeMap)

    #Change the subgraph to calculate H for the target node
    subgraph.σ[1] = stepInfo.targetCell
    subgraph.τ[1] = CPM.cell.types[stepInfo.targetCell]

    #Loop through all edges in network and calculate adhesion adjacencies
    targetH = AP(subgraph, CPM.cell.typeMap)

    #Calculate the adhesion difference between target and source
    return targetH - sourceH
end

function (VP::VolumePenalty)(CPM::CellPotts, stepInfo::MHStepInfo)

    #Create a vector of source and target memberships
    sourceTarget = [stepInfo.sourceCell, stepInfo.targetCell]

    #Extract cell types
    cellTypeIdx = [CPM.cell.typeMap[type] for type in CPM.cell.types[sourceTarget]] 

    #Extract the desired volumes
    desiredVolumes = VP.desiredVolumes[sourceTarget]

    #Extract the current volumes
    currentVolumes = CPM.cell.volumes[sourceTarget]

    #Calculate the target volumes
    targetVolumes = currentVolumes .+ [-1, 1]

    #Return the delta (Target - Source)
    return sum(@. VP.λᵥ[cellTypeIdx]*(targetVolumes - desiredVolumes)^2) - sum(@. VP.λᵥ[cellTypeIdx]*(currentVolumes - desiredVolumes)^2)
end


####################################################
# Helper Functions that need structs defined
####################################################

function UpdateConnections!(graph::NamedGraph; checkConnect::Bool=false)

    #Reset the articulation points (is this needed?)
    graph.isArticulation .= falses(size(graph.isArticulation))

    #Loop through cells to find articulation points
    for σᵢ in unique(graph.σ)
        
        #Get the subgraph for a given cell ID
        cellIdx = findall(isequal(σᵢ), graph.σ)
        subgraph = graph.network[cellIdx]

        if checkConnect
            if !is_connected(subgraph) #if not connected
                error("some cells are disconnected, try rerunning or use a different cell initialization")
            end
        end

        #Update the articulation points
        graph.isArticulation[cellIdx[articulation(subgraph)]] .= true
    end

    return nothing
end

#Overload getindex to get subset of Namedgraph
Base.getindex(graph::NamedGraph, i) = NamedGraph(graph.network[i], graph.σ[i], graph.τ[i], graph.isArticulation[i])

####################################################
# Functions to set initial cell locations
####################################################

function GrowCells(gr::SimpleGraph, M::ModelParameters)

    #initialize matrix of cell indicies (σ)
    cellMembership = zeros(Int, M.graphDimension...)
    
    #Find the center node
    centerIdx = CartesianIndex(M.graphDimension.÷2)
    nodeIdx = LinearIndices(M.graphDimension)
    centerNode = nodeIdx[centerIdx]

    #Determine how far nodes are from the center
    nodeDis = gdistances(gr, centerNode)

    #How many nodes need to be initialized?
    totalNodes = M.cellCounts ⋅ M.cellVolumes #dot product to sum and multiply

    #Get a sorted permutation of the distance
    sortedDis = sortperm(nodeDis)

    #Assign 1:totalNodes to be filled with cells and the rest medium
    networkIdx = sortedDis[1:totalNodes] 

    #Partition the identified nodes by the number of cells needed
    if sum(M.cellCounts) == 1 #There is only one cell (no need to partition)
        cellMembership[networkIdx] .= 1
    else
        cellMembership[networkIdx] = Metis.partition(gr[networkIdx], sum(M.cellCounts))
    end

    return cellMembership
end

####################################################
# Override Base.show for each struct
####################################################

function Base.show(io::IO, CPM::CellPotts) 
    println("Cell Potts Model:")
    #Grid
    dim = length(CPM.M.graphDimension)
    if dim == 1
        println("$(CPM.M.graphDimension[1])×$(CPM.M.graphDimension[1]) Grid")
    elseif dim == 2
        println("$(CPM.M.graphDimension[1])×$(CPM.M.graphDimension[2]) Grid")
    else
        println("$(CPM.M.graphDimension[1])×$(CPM.M.graphDimension[2])×$(CPM.M.graphDimension[3]) Grid")
    end

    #Cells and types
    print("Cell Counts:")
    allCellTypes = 
    for (key, value) in countmap(CPM.cell.types) #remove medium
        if key ≠ "Medium"
        print(" [$(key) → $(value)]")
        end
    end

    if length(CPM.M.cellTypes) > 1
        println(" [Total → $(length(CPM.cell.ids))]")
    else
        print("\n")
    end

    print("Model Penalties:")
    for p in typeof.(CPM.M.penalties)
        print(" $(p)")
    end
    print("\n")
    println("Current Energy: ", CPM.energy)
    println("Grid Temperature: ", CPM.M.temperature)
    println("Steps: ", CPM.stepCounter)
end


####################################################
# Metropolis–Hasting Step
####################################################

function MHStep!(CPM::CellPotts)

    #Create a structure to hold the source and target information
    stepInfo = MHStepInfo() #should add as input, might save a little time

    #Loop through until a good source target is found
    searching = true
    while searching
        #Pick a random location on the graph
        stepInfo.sourceNode = rand(1:nv(CPM.graph.network))
        #What cell does it belong to?
        stepInfo.sourceCell = CPM.graph.σ[stepInfo.sourceNode]

        #Get all of the unique cell IDs neighboring this Node
        stepInfo.sourceNodeNeighbors = neighbors(CPM.graph.network, stepInfo.sourceNode)
        possibleCellTargets = unique(CPM.graph.σ[stepInfo.sourceNodeNeighbors])
        #Choose a target
        stepInfo.targetCell = rand(possibleCellTargets)

        #Some checks before attempting a flip
            #In the middle of a cell
            inMiddle = all(possibleCellTargets .== stepInfo.sourceCell)
            #will fragment the cell
            isArticulation = CPM.graph.isArticulation[stepInfo.sourceNode]
            #target is the same as source cell 
            isSource = stepInfo.targetCell == stepInfo.sourceCell

        #if all checks pass, attempt flip
        if !(inMiddle | isArticulation | isSource) 
            searching = false
        else
            CPM.stepCounter += 1
        end
    end    

    #Calculate the change in energy when source node is modified
    ΔH =  sum([f(CPM, stepInfo) for f in CPM.M.penalties])

    #Calculate an acceptance ratio
    acceptRatio = min(1.0,exp(-ΔH/CPM.M.temperature))

    if rand() < acceptRatio #If the acceptance ratio is large enough

        #Need to update all cell and graph properties
        #---Cell properties---

        #Update cell volumes
        CPM.cell.volumes[stepInfo.sourceCell] -= 1
        CPM.cell.volumes[stepInfo.targetCell] += 1


        #---Graph properties---

        #Cell IDs
        CPM.graph.σ[stepInfo.sourceNode] = stepInfo.targetCell
        CPM.graph.τ[stepInfo.sourceNode] = CPM.cell.types[stepInfo.targetCell]

        #Articulation points
        UpdateConnections!(CPM.graph)
        
        #---Overall properties---
        #Update energy
        CPM.energy += ΔH
        #Update visual
        CPM.visual[stepInfo.sourceNode] = stepInfo.targetCell

    end
    return nothing
end

####################################################
# Cell Division
####################################################

#This is so easy with a graph 
function Division(CPM::CellPotts, σ::Int)

    cellIdx = findall(isequal(σ), CPM.graph.σ)

    #Vector of ones and twos
    newCellsIdx = Metis.partition(CPM.graph.network[cellIdx], 2)
    newCellNodeIdx = cellIdx[newCellsIdx .== 2]

    #Update cell attributes
    #ID
        push!(CPM.cell.ids, CPM.cell.ids[end]+1)
    #Volumes
        CPM.cell.volumes[σ] = sum(x->x==1, newCellsIdx)
        push!(CPM.cell.volumes, sum(x->x==2, newCellsIdx))
    #types
        push!(CPM.cell.types, CPM.cell.types[σ])

    #Update graph attributes
    #σ
        CPM.graph.σ[newCellNodeIdx] .= CPM.cell.ids[end]
    #τ
        CPM.graph.τ[newCellNodeIdx] .= CPM.cell.types[end]
    #isArticulation
        UpdateConnections!(CPM.graph)

    #Penalties???   

    #Global CellPotts attributes
    #energy
        CPM.energy = sum([f(CPM) for f in CPM.M.penalties])
    #visual
        CPM.visual[newCellNodeIdx] .= CPM.cell.ids[end]
end

####################################################
# Testing initialization
####################################################

# M = ModelParameters(
#     graphDimension = (100,500),
#     isPeriodic = true,
#     cellTypes = ["TCell", "BCell"],
#     cellCounts = [3, 5],
#     cellVolumes = [10, 8],
#     penalties = [AdhesionPenalty([0 1 1; 1 1 1; 1 1 1]),
#                  VolumePenalty(fill(9,18),[1,1])]
# )

# CPM = CellPotts(M)


M = ModelParameters(
    graphDimension = (50,50),
    isPeriodic = true,
    cellTypes = ["TCell"],
    cellCounts = [20],
    cellVolumes = [50],
    penalties = [AdhesionPenalty([0 50; 50 100]),
                 VolumePenalty(fill(50,20),[5])],
    temperature = 20.0)

CPM = CellPotts(M)

# using Plots
# gr()

# heatmap(CPM.visual,c = :Purples)


# for i in 1:100_000
#     MHStep!(CPM)
#     if mod(i,100) == 0
#         display(heatmap(CPM.visual, c = :Purples))
#         @show i
#     end
# end

####################################################
# Testing 3D
####################################################

# M = ModelParameters(
#     graphDimension = (40,40,40),
#     isPeriodic = true,
#     cellTypes = ["TCell"],
#     cellCounts = [20],
#     cellVolumes = [1000],
#     penalties = [AdhesionPenalty([0 50; 50 100]),
#                  VolumePenalty(fill(1000,20),[5])],
#     temperature = 20.0)
# CPM = CellPotts(M)
