abstract type Penalty end

#Offset arrays so the zeroth index refers to Medium

mutable struct AdhesionPenalty <: Penalty
    J::OffsetMatrix{Int, Matrix{Int}} #J[n,m] gives the adhesion penality for cells with types n and m

    function AdhesionPenalty(J::Matrix{Int})
        issymmetric(J) ? nothing : error("J needs to be symmetric")
        
        Joff = OffsetArray(J, CartesianIndex(0, 0):CartesianIndex(size(J).-1))
        return new(Joff)
    end
end

mutable struct VolumePenalty <: Penalty
    λᵥ::OffsetVector{Int,Vector{Int}}

    function VolumePenalty(λᵥ::Vector{Int})
        λᵥOff = OffsetVector([0; λᵥ], 0:length(λᵥ))
        return new(λᵥOff)
    end
end

mutable struct PerimeterPenalty <: Penalty
    λₚ::OffsetVector{Int,Vector{Int}}

    function PerimeterPenalty(λᵥ::Vector{Int}) 
        λᵥOff = OffsetVector([0; λᵥ], 0:length(λᵥ))

        return new(λᵥOff, currentPerimeters)
    end
end

#TODO Migration should work without floats
mutable struct MigrationPenalty{N} <: Penalty
    memory::Array{Int,N}
    maxAct::Int
    λ::Int

    function MigrationPenalty(maxAct::Int, λ::Int, gridSize::NTuple{N, Int}) where N
        memory = zeros(Int, gridSize)

        return new{N}(memory, maxAct, λ)
    end
end