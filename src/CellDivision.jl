
#Divide cell using slope and line technique
function CellDivide!(CPM::CellPotts, σ::Int64)
    
    #Check if cell crosses the grid boundary
    doesCross = (σ in CPM.grid[[1,CPM.n],:]) | (σ in CPM.grid[:,[1,CPM.n]])

    #Shift the grid so the cell isn't segmented
    if doesCross
        shiftGrid = gridShifter_lazy(CPM,σ)

        #Extract the cell from the grid (like zoom into it)
        corners = extrema(findall( shiftGrid .== σ))
        window = @view shiftGrid[ corners[1]:corners[2] ] #this is cool

    else #Extract normally
        corners = extrema(findall( CPM.grid .== σ))
        window = @view CPM.grid[ corners[1]:corners[2] ]
    
    end

    #Create a function that takes in a slope and outputs difference in daughter cell sizes
    DaughterSizeDiff = slope -> abs( diff(OptimizeSlope(slope, window, σ))[1] ) #L2 norm?
    opt = optimize(DaughterSizeDiff, -100.0, 100.0) #pick a slope b/w -100 and 100

    #store the volumes for the daughter cells form the best slope
    aboveBelow =  OptimizeSlope(opt.minimizer, window, σ)

    #Use the optimized slope to update the window
    CPM.σ += 1 #adding a new cell id
    LineDivider!(window, opt.minimizer, σ, CPM.σ) #this updates the grid

    ### Update the rest of CPM ###

    #volume desired for new cell
    push!(CPM.Vd, CPM.Vd[σ]) 
    #volume current for both cells
    push!(CPM.Vc, aboveBelow[1]) #new cell is above
    CPM.Vc[σ] = aboveBelow[2]
    #volume lagrange multiplier
    push!(CPM.λᵥ, CPM.λᵥ[σ]) 
    #Hamiltonian
    #make this a function?
        #H (adhesion)
        H = 0
        for I in CartesianIndices(CPM.grid)
            H += sum(Neighbors(I,CPM.n) .≠ CPM.grid[I]) #add 1 for each dis-similar neighbor
        end

        #H Volume
        H += sum(@. CPM.λᵥ*(CPM.Vc - CPM.Vd)^2) #difference b/w desired and current volume
    CPM.H = H

    return nothing
end

#Draw a line through the cell and count the number of square above and below
function OptimizeSlope(slope, window, σ)
    #Take a point in the middle of the window
    (midX,midY) = size(window).÷2
    #Use the slope to calculate the y-intercept
    b = midY - slope*midX

    #Where are the cell squares of interest
    posFinds = findall(window .== σ) #this repeat find could be avoided

    #Loop through the poisition and cound cell squares above/below line
    aboveBelow = [0, 0]
    for coord in posFinds
        if slope*coord[1] + b < coord[2]
            aboveBelow[1] +=1
        else
            aboveBelow[2] +=1
        end
    end

    #return the counts for above and below to minimize
    return aboveBelow 
end



#Draw a line through the cell and count the number of square above and below
function LineDivider!(window, slope, σ, σ_new)
    #Take a point in the middle of the window
    (midX,midY) = size(window).÷2
    #Use the slope to calculate the y-intercept
    b = midY - slope*midX

    #Where are the cell squares of interest
    posFinds = findall(window .== σ) #this repeat find could be avoided

    #Loop through the poisition and cound cell squares above/below line
    above = 0
    below = 0
    for coord in posFinds
        if slope*coord[1] + b < coord[2]
            window[coord] = σ_new
        end
    end
end

#Incrementally shift the array down and to the right until the cell is connected
function gridShifter_lazy(CPM::CellPotts,σ::Int64)
    #Amount of shifting to the left and right
    shifter = [0,0]

    #Create a view of the grid to shift around
    shiftGrid = @view CPM.grid[:,:]

    #Shift down until the cell is not on the borders
    while σ in shiftGrid[[1,CPM.n],:]
        shiftGrid = circle_shift(shiftGrid,[1,0])
        shifter[1] += 1
    end

    #Shift right until the cell is not on the borders
    while σ in shiftGrid[:,[1,CPM.n]]
        shiftGrid = circle_shift(shiftGrid,[0,1])
        shifter[2] += 1
    end

    #Return a lazily shifted grid
    return circle_shift(CPM.grid,shifter)
end