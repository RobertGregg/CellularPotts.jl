#This file contains functions that are pervasive throughout this package


#The mod1 function always returns a number between 1 and n
#0 gets mapped to n
#n+1 gets mapped to 1
#This is exactly what is needed for periodic boundaries
#Make mod1 work for CartesianIndex (I is the tuple of indices)
Base.mod1(x::CartesianIndex{2},y::Int)  = CartesianIndex(mod1.(x.I,y))


#Get the indices of the 8 (or 4) neighboring values
function Neighbors(loc::CartesianIndex{2}, n::Int64 ;type=8)
    Δx = CartesianIndex(1,0)
    Δy = CartesianIndex(0,1)

    if type==8
        return mod1.([loc-Δx-Δy,loc-Δy,loc+Δx-Δy,
                      loc-Δx,           loc+Δx,
                      loc-Δx+Δy,loc+Δy,loc+Δx+Δy],n)
    elseif type==4
        return mod1.([          loc-Δy,
                      loc-Δx,           loc+Δx,
                                loc+Δy         ],n)
    else
        throw(DomainError(type, "type must be 8 or 4"))
    end
end


#Create a structure for patrol movement
Base.@kwdef mutable struct PatrolMovement # specified default
    on::Bool = false #Do I want this in the simulation?
    Ma::Int64 = 20 #Maximum value for Gm
    λact::Float64 = 1.0 #act model lagrange multiplier
    Gm::Array{Int64,2} #Array to keep track of active cell patrol movement
end

#Create a structure to hold the model parameters with default parameters
mutable struct CellPotts
    n::Int64 #side length of the grid
    grid::Array{Int64,2} #array of cell ids denoting where cells are located
    β::Float64 #Simulation inverse temperature
    σ::Int64 #Number of unique IDs
    Vd::OffsetVector{Int64,Array{Int64,1}} #desired cell volumes
    Vc::OffsetVector{Int64,Array{Int64,1}} #current cell volumes
    λᵥ::OffsetVector{Int64,Array{Int64,1}} #volume lagrange multiplier
    H::Real #Hamiltonian energy (volume, adhesion)
    pm::PatrolMovement

    #= Note about Offset Arrays
        The medium (grid square with no cell) is given a value of zero in the grid.
        Sometimes a grid will attempt to change from a cell to medium (e.g. 4 → 0).
        It is convenient to use the grid id as an index (grid id of 1 corresponds with the volume vector's 1st index).
        This means I need an index of 0 to refer to the medium (hence the base 0 offset array)
    =#

    #Inner constructor
    #Loop through each grid point and count the number of different pixels
    #calculate the initial energy
    function CellPotts(;n=100,β=3.0,σ=20,Vd=rand(40:55,20),patrol=false)

        #Initialize the grid
        grid = rand(0:σ,n,n)

        #lagrange multipliers
        λᵥ = OffsetVector([0;ones(Int64,σ)],0:σ)

        #Vd (desired volumes)
        Vd = OffsetVector([0;Vd],0:σ) #the medium get index 0 and cell 1 gets index 1

        #Vc (initial current volume)
        Vc = OffsetVector(counts(grid),0:σ) #(0 means no cell on gridpoint)

        #Active cell memory grid (equal to Ma if grid square contains cell)
        pm = PatrolMovement(Gm = zeros(size(grid)))

        if patrol
            pm.on = patrol #change to true
            pm.Gm = @. pm.Ma*(grid ≠ 0) # if there is a cell id, replace 0 with Ma in Gm
        end

        #H (adhesion)
        H = 0
        for I in CartesianIndices(grid)
            H += sum(Neighbors(I,n) .≠ grid[I]) #add 1 for each dis-similar neighbor
        end

        #H Volume
        H += sum(@. λᵥ*(Vc - Vd)^2) #difference b/w desired and current volume


        #Return a new instantiation
        return new(n,grid,β,σ,Vd,Vc,λᵥ,H,pm)
    end
end


#Make the CellPotts struct print nicely
function Base.show(io::IO, c::CellPotts) 
    println("Cell Potts Model:")
    @printf "%d×%d grid with %d cells\n" c.n c.n c.σ
    println("Current Energy: ",c.H)
end


#Calculate the change in energy from adhesion
function Propose!(CPM::CellPotts,
                idCurrent::Int64,
                idPropose::Int64,
                neighborIndices::Array{CartesianIndex{2},1})

    #Choose a random neighbor ID
    neighborIDs = CPM.grid[neighborIndices]

    #Calculate ΔH from adhesion
    ΔH = sum(neighborIDs .≠ idPropose) - sum(neighborIDs .≠ idCurrent)

    Vpropose = copy(CPM.Vc) #new proposed volumes
    Vpropose[idCurrent] -= 1 #lower the current volume id
    Vpropose[idPropose] += 1 #increase the porposed volume id

    for (λ,Vd,Vc,Vp) in zip(CPM.λᵥ, CPM.Vd,CPM.Vc,Vpropose)
        ΔH += λ*((Vp - Vd)^2 - (Vc - Vd)^2) #This might be wrong
    end

    return ΔH
end

#Do a metropolis hastings step
function MHStep!(CPM::CellPotts)

    #Pick a random location and pick a neighboring location to copy from
    #make rand output tuple? e.g. CartesianIndex(Tuple(rand(1:10,2))...)
    #The tuple method is slower and causes allocations
    locCurrent = CartesianIndex(rand(1:CPM.n),rand(1:CPM.n))
    idCurrent =  CPM.grid[locCurrent]

    #Determine the neighbors of that random location
    neighborIndices = Neighbors(locCurrent,CPM.n)

    locPropose = rand(neighborIndices) #pick and random neighbor to get a new id
    idPropose =  CPM.grid[locPropose]

    ΔH = Propose!(CPM,idCurrent,idPropose,neighborIndices) #adhesive and volume changes

    if CPM.pm.on
        ΔH_Gm = GmStep!(CPM,locCurrent,locPropose)
        ΔH += ΔH_Gm
    end

    acceptRatio = min(1.0,exp(-ΔH*CPM.β))

    if (rand()<acceptRatio) & (idCurrent ≠ idPropose) #If we like the move update grid, else do nothing
        #Update the current volume
        CPM.Vc[idCurrent] -= 1 #lower the current volume id
        CPM.Vc[idPropose] += 1 #increase the porposed volume id

        #Update the grid
        CPM.grid[locCurrent] = idPropose

        #Update Gm
        if CPM.pm.on
            #Update the random proposal if it's part of a cell
            CPM.pm.Gm[locCurrent] = idPropose == 0 ? 0 : CPM.pm.Ma

            #Subtract 1 from all of Gm, except for 0 values
            for I in CartesianIndices(CPM.pm.Gm)
                if CPM.pm.Gm[I] > 0 
                    CPM.pm.Gm[I] -= 1
                end
            end
        end

        #Update H
        CPM.H += ΔH
    end

    return nothing
end

#Based off of this paper: https://doi.org/10.1371/journal.pcbi.1004280

function GmPropose!(CPM::CellPotts,loc::CartesianIndex{2})
    
    neighborIndices = push!(Neighbors(loc,CPM.n),loc) #take location and append neighbors

    #remove neighbors that do not have the same id
    filter!(x -> CPM.grid[x] == CPM.grid[loc],neighborIndices)

    return geomean(CPM.pm.Gm[neighborIndices])
end

function GmStep!(CPM::CellPotts,locCurrent::CartesianIndex{2},locPropose::CartesianIndex{2})

    #Now we have a location and a new location to try and extend to
    ΔH_Gm = (CPM.pm.λact / CPM.pm.Ma) * (GmPropose!(CPM,locPropose) - GmPropose!(CPM,locCurrent))
    
    return ΔH_Gm
end

#This is very ugly and maybe one day I'll make it better
#This function takes in the length of a square grid and returns pairs of all adjacent squares (wraps around)
# ╔═══╤═══╗
# ║ 1 │ 3 ║
# ╠═══╪═══╣
# ║ 2 │ 4 ║
# ╚═══╧═══╝
# Edge2Grid(2) = [[2, 1], [4, 3], [1, 2], [3, 4], [2, 1], [4, 3], [2, 4], [1, 3], [4, 2], [3, 1], [2, 4], [1, 3]]

function Edge2Grid(gridSize::Real)
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