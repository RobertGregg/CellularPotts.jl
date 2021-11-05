abstract type Home end

mutable struct Condo <: Home
    color::String
    price::Float64
end

struct Apartment <: Home
    color::String
    floors::Int64
    price::Float64
end


mutable struct Block
    homes::Vector{Home}
    totalPriceApartments::Float64
    totalPriceCondos::Float64
    totalPrice::Float64
end

myHomes = vcat([Condo("Red", 100_000.00*rand()) for _ in 1:10000], [Apartment("Blue", rand(1:3), 80_000.50*rand()) for _ in 1:10000])
myBlock = Block(myHomes, 0.0, 0.0, 0.0)



function calcPrice!(block::Block, home::Condo)
    block.totalPriceCondos += home.price
    return nothing
end

function calcPrice!(block::Block, home::Apartment)
    block.totalPriceApartments += home.price * home.floors #assume each floor has one apartment, just making this up
    return nothing
end



function setTotalPrice(A::Block)
    
    A.totalPriceApartments = A.totalPriceCondos = A.totalPrice = 0.0

    @inbounds for i in eachindex(A.homes)
        #this function barrier technique avoids slow-down from type instability and will not allocate if it returns nothing 
        calcPrice!(A, A.homes[i]) 
    end

    A.totalPrice = A.totalPriceApartments + A.totalPriceCondos

    return nothing
end



using BenchmarkTools

@benchmark setTotalPrice(myBlock)