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

        return new(λᵥOff)
    end
end

#TODO How to apply penalties to select cells?
struct MigrationPenalty <: Penalty
    maxAct::Int
    λ::Int
    cellTypes::Vector{Symbol}
end