#Custom sparsevector that does not allocate

mutable struct SV{T<:Integer} <: AbstractVector{T}
    n::T
    nzind::Vector{T}
    nzval::Vector{T}

    function SV(n::T) where {T<:Integer}
        nzind = T[]
        sizehint!(nzind,n)

        nzval = T[]
        sizehint!(nzval,n)
        return new{T}(n,nzind,nzval)
    end
end

Base.size(S::SV) = (S.n,)
Base.IndexStyle(::Type{<:SV}) = IndexLinear()

nnz(S::SV) = length(S.nzind)

function Base.getindex(S::SV, i::Int)
    ii = searchsortedfirst(S.nzind,i)
    
    return (ii ≤ nnz(S) && S.nzind[ii] == i) ? S.nzval[ii] : 0
end


function Base.setindex!(S::SV{T}, val::T, i::T) where T <: Integer
    checkbounds(S, i)

    if iszero(val)
        return nothing
    end

    ii = searchsortedfirst(S.nzind,i)

    if ii ≤ nnz(S) && S.nzind[ii] == i
        S.nzval[ii] = val
    else
        insert!(S.nzind, ii, i)
        insert!(S.nzval, ii, val)
    end
end


s = SV(10)

s.nzind = [1,4]
s.nzval = [10,20]

s[3] = 100

sum(s)



dict = Dict{Int, Int}()
sizehint!(dict, 200)
@benchmark begin
    for i in 1:100
        dict[i] = i^2
    end
end