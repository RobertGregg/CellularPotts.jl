using Plots
using CellularPotts

#TODO periodic boundaries, diagonals not getting catched 

#For diagonals you can start on a bit sticking out and it will just loop around and halt
#maybe do something like push (NaN,NaN) to  prevent the while stop condition?

#=
Iterates down columns first then over to the right for rows
#Number the cases by their binary representation
+-----------|----------+        
|           |          |        
|    0      |    0     |        
|         +---+        |        
----------|---|--------- ------> [0,1,0,1] --> case 5 (top straight edge)
|         +---+        |        
|    1      |    1     |        
|           |          |        
+-----------|----------+      

We're always moving counter-clockwise, so for case 5 we go to the left

#Corners are referenced by the bottom right node (the 1's place in binary vector)
=#

Base.mod1(I::CartesianIndex{N}, i::NTuple{N, Int64}) where N = CartesianIndex(mod1.(Tuple(I),i))
getCorner(I) = Tuple(I) .- 0.5
checkNode(M,cell,I, isPeriodic) = (isPeriodic || I ∈ CartesianIndices(M)) && isequal(M[I],cell)

function createShape(cpm,cell)

    #get the matrix of IDs
    #transpose needed b/c of how plots displays matrices
    M = cpm.space.nodeIDs'

    #Is the space periodic?
    isPeriodic = cpm.space.periodic

    #These could be const or put into some struct
    #Save a vector where indexes are cases and entries are directions to go
    up = CartesianIndex(-1,0)
    left = CartesianIndex(0,-1)
    down = CartesianIndex(1,0)
    right = CartesianIndex(0,1)

    #Looking like the konami code
    direction = [
        down,  #1  [0 0; 0 1] bottom right corner 
        right, #2  [0 1; 0 0] top right corner
        down,  #3  [0 1; 0 1] right edge
        left,  #4  [0 0; 1 0] bottom left corner
        left,  #5  [0 0; 1 1] bottom edge
        left,  #6  [0 1; 1 0] diagonal flipped
        left,  #7  [0 1; 1 1] bottom right wedge
        up,    #8  [1 0; 0 0] top left corner
        up,    #9  [1 0; 0 1] diagonal
        right, #10 [1 1; 0 0] top edge
        down,  #11 [1 1; 0 1] top right wedge
        up,    #12 [1 0; 1 0] left edge
        up,    #13 [1 0; 1 1] bottom left wedge
        right  #14 [1 1; 1 0] top left wedge
    ]
   
    corners = Tuple{Float64,Float64}[]
    I = findfirst(isequal(cell),M)
    J = I
    shapeConnected = false

    while !shapeConnected
        #add corner to list
        push!(corners, getCorner(J))
    
        #Go to the next corner by calculating which case you have
        #Need to check the edges
        # case = M[J] + 2*M[J + up] + 4*M[J + left] + 8*M[J + up + left]
        case = 0
    
        if checkNode(M, cell, J, isPeriodic) 
            case += 1
        end
    
        if checkNode(M, cell, J + up, isPeriodic) 
            case += 2
        end
    
        if checkNode(M, cell, J + left, isPeriodic) 
            case += 4
        end
    
        if checkNode(M, cell, J + up + left, isPeriodic) 
            case += 8
        end
            
        J += direction[case]
    
        # if J ∉ CartesianIndices(M)
        #     push!(corners,(NaN,NaN))
        # end

        J = mod1(J, size(M))

        if I == J
            push!(corners,getCorner(I))
            shapeConnected = true
        end
    end

    return corners
end


#############################################
# testing
#############################################

cpm = CellPotts(
        CellSpace(10,10, periodic=false, diagonal=false),
        CellState(names=[:A], volumes=[10], counts=[1]),
        [AdhesionPenalty(fill(10,2,2)), VolumePenalty([2])]
    )

ModelStep!(cpm)
visualize(cpm)

corners = [createShape(cpm,i) for i in 1:countcells(cpm)]

s = Shape.(corners)

heatmap(cpm.space.nodeIDs, c = :berlin)
plot!(s; fill=[(r,1.0) for r in cgrad(:rainbow1,countcells(cpm),categorical=true)])


scatter!(corners, markersize=1.5)
annotate!([(coord..., text(string(coord),4)) for coord in filter(!isequal((NaN,NaN)), corners[1])])
annotate!([(coord..., text(string(coord),4)) for coord in corners])

 

