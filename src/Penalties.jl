abstract type Penalty end

#Offset arrays so the zeroth index refers to Medium

struct AdhesionPenalty <: Penalty
    J::OffsetMatrix{Int, Matrix{Int}} #J[n,m] gives the adhesion penality for cells with types n and m

    function AdhesionPenalty(J::Matrix{Int})
        issymmetric(J) ? nothing : error("J needs to be symmetric")
        
        Joff = OffsetArray(J, Origin(0))
        return new(Joff)
    end
end

struct VolumePenalty <: Penalty
    λᵥ::OffsetVector{Int,Vector{Int}}

    function VolumePenalty(λᵥ::Vector{Int})
        λᵥOff = OffsetVector([0; λᵥ], Origin(0))
        return new(λᵥOff)
    end
end

struct PerimeterPenalty <: Penalty
    λₚ::OffsetVector{Int,Vector{Int}}

    function PerimeterPenalty(λᵥ::Vector{Int}) 
        λᵥOff = OffsetVector([0; λᵥ], Origin(0))

        return new(λᵥOff, currentPerimeters)
    end
end

#TODO Migration should work without floats
struct MigrationPenalty{N} <: Penalty
    memory::Array{Int,N}
    maxAct::Int
    λ::Int

    function MigrationPenalty(maxAct::Int, λ::Int, gridSize::NTuple{N, Int}) where N
        memory = zeros(Int, gridSize)

        return new{N}(memory, maxAct, λ)
    end
end