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
    end
end


#Create a structure to hold the model parameters with default parameters
mutable struct CellPotts
    n::Int64 #side length of the grid
    grid::Array{Int64,2} #array of cell ids denoting where cells are located
    β::Float64 #Simulation inverse temperature
    σ::Int64 #Number of unique IDs
    Vd::Vector{Int64} #desired cell volumes
    Vc::Vector{Int64} #current cell volumes
    λ::Vector{Int64} #volume lagrange multiplier
    H::Int64 #Hamiltonian energy (volume, adhesion)

    #Inner constructor
    #Loop through each grid point and count the number of different pixels
    #calculate the initial energy
    function CellPotts(;n=100,β=3.0,σ=20,Vd=vcat(0,rand(40:55,19)))

        #Initialize the grid
        grid = rand(1:σ,n,n)

        #lagrange multipliers
        λ = vcat(0,ones(Int64,σ-1))

        #Vc (initial current volume)
        Vc = counts(grid) #first count (0 means no cell on gridpoint)

        #H (adhesion)
        H = 0
        for I in CartesianIndices(grid)
            H += sum(Neighbors(I,n) .≠ grid[I]) #add 1 for each dis-similar neighbor
        end

        #H Volume
        H += sum(@. λ*(Vc - Vd)^2) #difference b/w desired and current volume

        #Return a new instantiation
        return new(n,grid,β,σ,Vd,Vc,λ,H)
    end
end


#Make the CellPotts struct print nicely
function Base.show(io::IO, c::CellPotts) 
    println("Cell Potts Model:")
    @printf "%d×%d grid with %d unique cell types\n" c.n c.n c.σ
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
        ΔH += λ*((Vp - Vd)^2 - (Vc - Vd)^2)
    end

    return ΔH
end

#Do a metropolis hastings step
function MHStep!(CPM::CellPotts)

    #Pick a random number and put it in a random spot
    #make rand output tuple? e.g. CartesianIndex(Tuple(rand(1:10,2))...)
    #The tuple method is slower and causes allocations
    loc = CartesianIndex(rand(1:CPM.n),rand(1:CPM.n))
    id = rand(1:CPM.σ)

    ΔH = Propose!(CPM,loc,id)

    acceptRatio = min(1,exp(-ΔH*CPM.β))

    if rand()<acceptRatio #If we like the move update grid, else do nothing
        #Update the current volume
        CPM.Vc[CPM.grid[loc]] -= 1 #lower the current volume id
        CPM.Vc[id] += 1 #increase the porposed volume id

        #Update the grid
        CPM.grid[loc] = id
        #Update H
        CPM.H += ΔH
    end

    return nothing
end