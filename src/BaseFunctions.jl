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


#Create a structure to hold the model parameters with default parameters
mutable struct CellPotts
    n::Int64 #side length of the grid
    grid::Array{Int64,2} #array of cell ids denoting where cells are located
    β::Float64 #Simulation inverse temperature
    σ::Int64 #Number of unique IDs
    Vd::OffsetVector{Int64,Array{Int64,1}} #desired cell volumes
    Vc::OffsetVector{Int64,Array{Int64,1}} #current cell volumes
    λ::OffsetVector{Int64,Array{Int64,1}} #volume lagrange multiplier
    H::Real #Hamiltonian energy (volume, adhesion)
    Gm::Array{Int64,2} #array to keep track of active cell patrol movement
    Ma::Int64 #Maximum value for Gm

    #= Note about Offset Arrays
        The medium (grid square with no cell) is given a value of zero in the grid.
        Sometimes a grid will attempt to change from a cell to medium (e.g. 4 → 0).
        It is convenient to use the grid id as an index (grid id of 1 corresponds with the volume vector's 1st index).
        This means I need an index of 0 to refer to the medium (hence the base 0 offset array)
    =#

    #Inner constructor
    #Loop through each grid point and count the number of different pixels
    #calculate the initial energy
    function CellPotts(;n=100,β=3.0,σ=20,Vd=rand(40:55,20),patrol=true)

        #Initialize the grid
        grid = rand(0:σ,n,n)

        #lagrange multipliers
        λ = OffsetVector([0;ones(Int64,σ)],0:σ)

        #Vd (desired volumes)
        Vd = OffsetVector([0;Vd],0:σ) #the medium get index 0 and cell 1 gets index 1

        #Vc (initial current volume)
        Vc = OffsetVector(counts(grid),0:σ) #(0 means no cell on gridpoint)

        #Active cell memory grid (equal to Ma if grid square contains cell)
        if patrol
            Ma = 20
            Gm = @. Ma*(grid ≠ 0)
        else
            Ma = 0
            Gm = zeros(Int64,size(grid))
        end

        #H (adhesion)
        H = 0
        for I in CartesianIndices(grid)
            H += sum(Neighbors(I,n) .≠ grid[I]) #add 1 for each dis-similar neighbor
        end

        #H Volume
        H += sum(@. λ*(Vc - Vd)^2) #difference b/w desired and current volume


        #Return a new instantiation
        return new(n,grid,β,σ,Vd,Vc,λ,H,Gm,Ma)
    end
end


#Make the CellPotts struct print nicely
function Base.show(io::IO, c::CellPotts) 
    println("Cell Potts Model:")
    @printf "%d×%d grid with %d cells\n" c.n c.n c.σ
    println("Current Energy: ",c.H)
end


#Calculate the change in energy from adhesion
function Propose!(CPM::CellPotts,loc::CartesianIndex{2},id::Int64)

    #Calculate ΔH from adhesion
    adjacent = CPM.grid[Neighbors(loc,CPM.n)]
    ΔH = sum(adjacent .≠ id) - sum(adjacent .≠ CPM.grid[loc])


    Vpropose = copy(CPM.Vc) #new proposed volumes
    Vpropose[CPM.grid[loc]] -= 1 #lower the current volume id
    Vpropose[id] += 1 #increase the porposed volume id

    for (λ,Vd,Vc,Vp) in zip(CPM.λ, CPM.Vd,CPM.Vc,Vpropose)
        ΔH += λ*((Vp - Vd)^2 - (Vc - Vd)^2) #This might be wrong
    end

    return ΔH
end

#Do a metropolis hastings step
function MHStep!(CPM::CellPotts)

    #Pick a random number and put it in a random spot
    #make rand output tuple? e.g. CartesianIndex(Tuple(rand(1:10,2))...)
    #The tuple method is slower and causes allocations
    loc = CartesianIndex(rand(1:CPM.n),rand(1:CPM.n))
    id = rand(0:CPM.σ)

    ΔH = Propose!(CPM,loc,id)

    if CPM.Ma ≠ 0 #maybe should have a better check for this
        ΔH_Gm, loc_gm = GmStep!(CPM)
        ΔH += ΔH_Gm
    end

    acceptRatio = min(1,exp(-ΔH*CPM.β))

    if rand()<acceptRatio #If we like the move update grid, else do nothing
        #Update the current volume
        CPM.Vc[CPM.grid[loc]] -= 1 #lower the current volume id
        CPM.Vc[id] += 1 #increase the porposed volume id

        #Update the grid
        CPM.grid[loc] = id

        #Update Gm
        if CPM.Ma ≠ 0
            #Update the random proposal
            CPM.Gm[loc] = id == 0 ? 0 : CPM.Ma

            #Update the active movement
            CPM.Gm[loc_gm] = CPM.Ma

            #Subtract 1 from all of Gm, except for 0 values
            @. CPM.Gm -= 1
            @. CPM.Gm[CPM.Gm < 0] = 0 #replace negative numbers with 0
        end

        #Update H
        CPM.H += ΔH
    end

    return nothing
end

#Based off of this paper: https://doi.org/10.1371/journal.pcbi.1004280

function GmPropose!(CPM::CellPotts,loc::CartesianIndex{2})
    
    adjacent = push!(Neighbors(loc,CPM.n),loc) #take location and append neighbors

    #remove neighbors that do not have the same id
    filter!(x -> CPM.grid[x] == CPM.grid[loc],adjacent)

    return geomean(CPM.Gm[adjacent])
end

function GmStep!(CPM::CellPotts)

    #Check if random grid point is on the border (if not pick a new grid point)
    notSureIfBorder = true

    #Pick a random grid point
    #local variable issue (pre-define before loop) ?
    loc = CartesianIndex(1,1)
    proposedLoc = CartesianIndex(1,1)

    while notSureIfBorder
        
        #Pick a random location
        loc = CartesianIndex(rand(1:CPM.n),rand(1:CPM.n))

        #Find neighboring squares that are different
        neighborhood = Neighbors(loc,CPM.n)
        outsideCell = findall(x -> CPM.grid[x] ≠ CPM.grid[loc],neighborhood)

        #If there is a different square, mark it to try and extend cell
        if ~isempty(outsideCell)
            proposedLoc = rand(neighborhood[outsideCell])
            notSureIfBorder = false
        end
    end

    #Now we have a location and a new location to try and extend to
    ΔH_Gm = (1.0/CPM.Ma) * (GmPropose!(CPM,loc) - GmPropose!(CPM,proposedLoc)) #assume λact = 1 for Now
    
    return (ΔH_Gm,proposedLoc)
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
