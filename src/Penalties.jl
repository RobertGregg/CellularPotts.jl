abstract type Penalty end

#Offset arrays so the zeroth index referes to Medium

mutable struct AdhesionPenalty <: Penalty
    J::Matrix{Int}

    function AdhesionPenalty(J::Matrix{Int})
        issymmetric(J) && error("J needs to be symmetric")
        return new(J)
    end
end

mutable struct VolumePenalty <: Penalty
    λᵥ::Vector{Int}
end

mutable struct PerimeterPenalty <: Penalty
    λₚ::Vector{Int}
end