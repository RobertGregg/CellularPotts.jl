using Plots
using CellularPotts


#############################################
# Version 2
#############################################

#= Notes
 - When plotting a matrix, the y-axis is flipped b/c rows are counted down and the y axis goes up
 - check out yflip=true when plotting
 - Okay now matrix coordinates and 
 - Lets start with ignoring boundaries b/c they are hard
=#


M = [
    0 0 0 1 1;
    0 1 1 1 1;
    0 0 1 0 0;
    0 0 1 0 0;
    0 0 0 0 0
]


#Reverse b/c CartesianIndex(ycoord, xcoord) when plotting
getCorner(I) = reverse(Tuple(I) .- 0.5)
checkNode(M,cell,I) = I ∈ CartesianIndices(M) && isequal(M[I],cell)

function createShape(M,cell)

    #These could be const or put into some struct
    #Save a vector where indexes are cases and entries are directions to go
    up = CartesianIndex(-1,0)
    left = CartesianIndex(0,-1)
    down = CartesianIndex(1,0)
    right = CartesianIndex(0,1)

    #Looking like the konami code
    direction = Dict(
    (false, false, false, true) => down,  #1  bottom right corner 
    (false, true, false, false) => right, #2  top right corner
    (false, true, false, true) => down,   #3  right edge
    (false, false, true, false) => left,  #4  bottom left corner
    (false, false, true, true) => left,   #5  bottom edge
    (false, true, true, false) => left,   #6  diagonal flipped
    (false, true, true, true) => left,    #7  bottom right wedge
    (true, false, false, false) => up,    #8  top left corner
    (true, false, false, true) => up,     #9  diagonal
    (true, true, false, false) => right,  #10 top edge
    (true, true, false, true) => down,    #11 top right wedge
    (true, false, true, false) => up,     #12 left edge
    (true, false, true, true) => up,      #13 bottom left wedge
    (true, true, true, false) => right    #14 top left wedge
    )
   
    corners = Tuple{Float64,Float64}[]
    I = findfirst(isequal(cell),M)
    J = I
    shapeConnected = false

    while !shapeConnected
        #add corner to list
        push!(corners, getCorner(J))
    
        #Go to the next corner by calculating which case you have
        case = (
            checkNode(M, cell, J + up + left),
            checkNode(M, cell, J + up),
            checkNode(M, cell, J + left),
            checkNode(M, cell, J)
        )

        J += direction[case]
    
        if I == J
            push!(corners,getCorner(I))
            shapeConnected = true
        end
    end

    return corners
end

corners = createShape(M,1)

s = Shape(corners)

plot(;
    yflip=true,
    c=:berlin,
    size = (500,500),
    axis = nothing,
    framestyle=:none,
    xlims = (0.5, size(M,2)+0.5),
    ylims = (0.5, size(M,1)+0.5),
    aspect_ratio = :equal,
    legend = false)

    
    
plot!(s; fillcolor = :steelblue, linewidth = 1.5, yflip=true)
    
plot!(Shape([(0.5,0.5),(size(M,2)+0.5,0.5),(size(M,2)+0.5,size(M,1)+0.5),(0.5,size(M,1)+0.5)]),
            fillcolor = nothing,     # Fill color of the shape
            linecolor = :black,      # Border color of the shape
            linewidth = 4,           # Border thickness
            )

scatter!(corners, markersize=1.5)
annotate!([(coord..., text(string(coord),6)) for coord in corners])



#############################################
# Version 1
#############################################

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

 


x = [0, 1, 1, 0]
y = [0, 0, 1, 1]
myshape = Shape(x, y)

plot(myshape,
            fillcolor = :orange,     # Fill color of the shape
            linecolor = :blue,      # Border color of the shape
            linewidth = 3,           # Border thickness
                 legend = false,
                      aspect_ratio = :equal)