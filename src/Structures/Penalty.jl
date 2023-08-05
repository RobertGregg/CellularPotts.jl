####################################################
# Penalties
####################################################

"""
    Penalty
An abstract type representing a constraint imposed onto the cellular potts model.

To add a new penalty, a new struct subtyping `Penalty` needs to be defined and the `addPenalty!()` function needs to be extended to include the new penalty.

**Note**: variables associated with a new penalty may need to be offset such that index 0 maps to :Medium, index 1 maps to :Cell1, etc.
"""
abstract type Penalty end

"""
    AdhesionPenalty(J::Matrix{Int})
A concrete type that penalizes neighboring grid locations from different cells.

Requires a symmetric matrix `J` where `J[n,m]` gives the adhesion penality for cells with types n and m. `J` is zero-indexed meaning `J[0,1]` and `J[1,0]` corresponds to the `:Medium` ↔ `:Cell1` adhesion penalty.

**Note**: `J` is automatically transformed to be a zero-indexed offset array.
"""
struct AdhesionPenalty{T} <: Penalty
    J::T

    function AdhesionPenalty(J::AbstractMatrix{<:Real})
        issymmetric(J) ? nothing : error("J needs to be symmetric")
        Joff = offset(J)

        T = typeof(Joff)
        return new{T}(Joff)
    end
end

"""
    VolumePenalty(λᵥ::Vector{Int})
A concrete type that penalizes cells that deviate from their desired volume.

Requires a vector `λᵥ` with n penalties where n is the number of cell types. `λᵥ` is zero-indexed meaning `λᵥ[0]` corresponds to the `:Medium` volume penalty (which is set to zero).

**Note**: `λᵥ` is automatically transformed to be a zero-indexed offset array and does not require the volume penalty for `:Medium`.
"""
struct VolumePenalty{T} <: Penalty
    λᵥ::T

    function VolumePenalty(λᵥ::AbstractVector{<:Real})
        λᵥOff = offset([0; λᵥ])

        T = typeof(λᵥOff)
        return new{T}(λᵥOff)
    end
end

"""
    PerimeterPenalty(λᵥ::Vector{Int})
A concrete type that penalizes cells that deviate from their desired perimeter.

Requires a vector `λₚ` with n penalties where n is the number of cell types. `λₚ` is zero-indexed meaning `λₚ[0]` corresponds to the `:Medium` perimeter penalty (which is set to zero).

**Note**: `λₚ` is automatically transformed to be a zero-indexed offset array and does not require the perimeter penalty for `:Medium`.
"""
mutable struct PerimeterPenalty{T} <: Penalty
    λₚ::T
    Δpᵢ::Int
    Δpⱼ::Int

    function PerimeterPenalty(λₚ::AbstractVector{<:Real}) 
        λₚOff = offset([0; λₚ])

        T = typeof(λₚOff)
        return new{T}(λₚOff, 0, 0)
    end
end

"""
    MigrationPenalty(maxAct, λ, gridSize)
A concrete type that encourages cells to protude and drag themselves forward.

Two integer parameters control how cells protude:
 - `maxAct`: A maximum activity a grid location can have
 - `λ`: A parameter that controls the strength of this penalty
 - 'gridSize': The size of the space, simply supply size(space)

Increasing `maxAct` will cause grid locations to more likely protrude. Increasing `λ` will cause those protusions to reach farther away. 
"""
mutable struct MigrationPenalty{T} <: Penalty
    maxAct::Int
    λ::T
    nodeMemory::SparseMatrixCSC{Int,Int}

    function MigrationPenalty(maxAct::I, λ::AbstractVector{<:Real}, gridSize::NTuple{N,I}) where {I<:Integer, N}
        λOff = offset([0; λ])

        T=typeof(λOff)
        return new{T}(maxAct, λOff, spzeros(T,gridSize))
    end
end


"""
    ChemoTaxisPenalty(λ, Species)
A concrete type that encourages cells to move up or down a concentration gradient.

Two integer parameters control how cells protude:
 - `λ`: A parameter that controls the strength of this penalty
 - `Species`: The concentration profile for a species that should match the size of the cell space

Species concentration profile can be updated dynamically (e.g. by an ODE)

Supplying a positive λ will move cells up the gradient, negative values down the gradient.
"""
mutable struct ChemoTaxisPenalty{T,M} <: Penalty
    λ::T
    species::M

    function ChemoTaxisPenalty(λ::AbstractVector{<:Real}, species::M) where M<:AbstractArray
        λOff = offset([0; λ])

        T = typeof(λOff)
        return new{T,M}(λOff, species)
    end
end


####################################################
# Penalty names
####################################################

#Given a penalty type, output a string for that type
#These make printing info for CPM easier.
penaltyType(::AdhesionPenalty) = "Adhesion"
penaltyType(::VolumePenalty) = "Volume"
penaltyType(::PerimeterPenalty) = "Perimeter"
penaltyType(::MigrationPenalty) = "Migration"
penaltyType(::ChemoTaxisPenalty) = "ChemoTaxis"

####################################################
# Override Base.show
####################################################

show(io::IO, P::Penalty) = print(io, penaltyType(P))