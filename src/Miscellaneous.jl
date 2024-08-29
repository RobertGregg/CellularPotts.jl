####################################################
# Small helpful functions used throughout
####################################################

#Kronecker delta function (used in AdhesionPenalty)
δ(x, y) = isequal(x,y) ? 1 : 0


#Returns a zero indexed array
offset(x) = OffsetArray(x, OffsetArrays.Origin(0))


#see https://github.com/JuliaArrays/OffsetArrays.jl/pull/137
#This is enough for this package b/c we only use 0-indexed vectors
deleteat!(v::OffsetVector, i::Integer) = deleteat!(parent(v), i-first(v.offsets))


#loop through all pairs in v
allpairs(v) = Iterators.filter(i -> isless(i...), Iterators.product(v,v))


#decode(["a","b"], [2,3]) = ["a","a","b","b","b"]
#Could be faster with manual loop but not a critical step
decode(values,counts) = reduce(vcat,fill.(values,counts))


function encode(x::AbstractVector{T}) where T
    d = Dict{T,Int}()

    for i in x
        if haskey(d,i)
            d[i] += 1
        else
            push!(d, i => 1)
        end
    end

    return d
end

#Given a desired cell volume, calculate minimum perimeter on a square lattice
#There are 3 pages of notes behind this equation
#It's related to the minimum perimeter for a polyomino which is 2ceil(2√V)
#TODO Honestly why isn't perimeter the literal perimeter?
estPerimeter(V::Int) = iszero(V) ? 0 : 4ceil(Int,2sqrt(V)-3) + 2ceil(Int,2sqrt(V+1)-4) + 14

