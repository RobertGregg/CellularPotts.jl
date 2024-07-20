# What are Cellular Potts Models?

For the following explanation we will be using this Cellular Potts Model:

```@example CPMmodel
using CellularPotts, Graphs, Random

Random.seed!(314)

cpm = CellPotts(
    CellSpace(10,10; periodic=false, diagonal=false),
    CellState(:Epithelial, 20,3,positions=[(4,6),(8,8),(8,3)]),
    [AdhesionPenalty([0 30; 30 0]),VolumePenalty([5])]
    )
```

# Background

Cellular Potts Models (CPMs) are essentially just non-negative integer matrices with very specific rules and interpretations. A value of zero represents a location in space where no cell is present (which we call "Medium") whereas a nonzero value denotes part of a cell.

```@setup CPMmodel
using CellularPotts, Plots, Random

Random.seed!(314)

cpm = CellPotts(
    CellSpace(10,10;periodic=false),
    CellState(:Epithelial, 20,3,positions=[(4,6),(8,8),(8,3)]), [AdhesionPenalty([0 30; 30 0]),VolumePenalty([5])]
    )

    space = copy(cpm.space.nodeIDs')

(rows,columns) = size(cpm.space) # hide

function plotSpace(space; numbers=true, gridlines=true, kwargs...)

    plotObject = heatmap(
        space,
        c = cgrad(:tol_light, rev=true),
        axis=nothing,
        legend=nothing,
        framestyle=:box,
        aspect_ratio=:equal,
        xlims=(0.5, rows+0.5),
        ylims=(0.5, columns+0.5),
        clim=(0,maximum(space));
        kwargs...
        )

    if numbers
        annotate!(plotObject,repeat(1:rows,inner=columns),repeat(1:columns,rows),vec(space))
    end

    if gridlines
        vline!(plotObject,0.5:(rows+0.5), c=:black)
        hline!(plotObject,0.5:(columns+0.5), c=:black)
    end

    return plotObject
end

function addBox!(plotObject,x,y;boxsize=95,boxcolor=RGBA(1,0.0,0,0.3))
    verts = [(-0.5, -0.5),(-0.5, 0.5),(0.5, 0.5),(0.5, -0.5),(-0.5, -0.5)]
    plot!(plotObject,[x],[y], marker = (Shape(verts), boxsize,boxcolor))
end

```

```@example CPMmodel
plotSpace(space;size = (500,500)) # hide
```



Here we have a 10$\times$10 grid with 3 cells all labeled "Epithelial". Note that cells should be connected; for example, all grid entries with a value of 1 are adjacent to each other. In `CellularPotts.jl` we can choose whether diagonal grid sites are considered adjacent with the `diagonal` keyword in the [`CellSpace`](@ref) function.

# Step Model Forward

CPMs take discrete time steps to update the cells which is handled by the `MHStep!` function in `CellularPotts.jl`. The process of updating takes three basic steps:

1. A random grid location called the "source" and an adjacent grid site called "target" is chosen. 
2. The model attempts to update source to match the target.
3. The update is accepted or rejected, given a calculated acceptance rate. 

In the example below, the proposal is accepted.

```@example CPMmodel
p1 = plotSpace(space;numbers=false,size = (1000,300)) # hide
annotate!(p1,5,7,"S") # hide
title!(p1,"Choose a location") # hide

p2 = plotSpace(space;numbers=false) # hide
addBox!(p2,5,7;boxsize=34) # hide
annotate!(p2,5,7,"S") # hide
annotate!(p2,5.5,7,"←") # hide
annotate!(p2,6,7,"T") # hide
title!(p2,"Attempt to Update") # hide

p3 = plotSpace(space;numbers=false) # hide
addBox!(p3,5,7;boxsize=34,boxcolor=RGB(249/255,187/255,170/255)) # hide
annotate!(p3,5,7,"✓") # hide
title!(p3,"Accept or Reject") # hide

plot(p1,p2,p3,layout = (1, 3)) # hide
```

If the source and target are from the same cell, then no update will occur. However, the source or target can be "Medium" possibly leading to the cell increasing or decreasing in size. 

Typically, $n$ proposals are attempted where $n$ is the total number of grid sites (100 in this case). Once all proposals are accepted or rejected, the model has taken one step forward.

# Penalties

The acceptance probability is governed by an "energy function" $H$ which is continually trying to be minimized. In its simplest form, there are two terms or "penalties" that contribute to the energy function: the [AdhesionPenalty](@ref) and [VolumePenalty](@ref).

## Adhesion

The adhesion penalty encourages adjacent grid sites to have the same value and contributes a positive penalty if there is a mismatch. The strength of the penalty is determined by a symmetric matrix $J$ where the $i^{th}$ and $j^{th}$ entry give the contact penalty for the $i^{th}$ and $j^{th}$ **cell type**. In our example J might looks like this:

```math
J = \begin{bmatrix}
0 & 30\\
30 & 0
\end{bmatrix}
```

Here, `J[1,1]` represents the Medium-Medium penalty, `J[1,2]` and `J[2,1]` represents the Cell-Medium penalty, `J[2,2]` and represents the Cell-Cell penalty. Note that we are assuming our cells are all the same type (e.g., Epithelial). If Cell 3 was labeled a T-Cell, for example, then we would add another row and column to $J$. 

To calculate the total adhesion penalty, we sum over all neighbors of all grid sites. 

```@example CPMmodel
function caculateAdhesion(cpm)

    #Track total adhesion
    totalAdhesion = 0

    #Get the adhesion penalty object
    AP = cpm.penalties[1]

    #Matrix of cell IDs (0,1,2,3) and 
    cellIDs = cpm.space.nodeIDs
    #Matrix cell types (0 → Medium or 1 → Epithelial)
    cellTypes = cpm.space.nodeTypes

    #Loop through each grid site
    for i in LinearIndices(cellIDs)

        #What cell is grid site i?
        σᵢ = cellIDs[i]
        #What type is grid site i?
        τᵢ = cellTypes[i]

        #Loop through neighbors
        for j in neighbors(cpm.space,i)
            τⱼ = cellTypes[j]
            σⱼ = cellIDs[j]
            #Check if the grid sites are from different cells
            if σᵢ ≠ σⱼ
                totalAdhesion += AP.J[τᵢ,τⱼ]
            end
        end
    end

    return(totalAdhesion)
end

println("There are $(caculateAdhesion(cpm)÷30) adhesion penalties") 
```

In practice, we don't need to calculate the total adhesion every time because it only changes locally around the source and target grid sites. 

The full equation for adhesion would look like this:

```math
\text{Total Adhesion Penalty} = \sum_{\substack{i,j \text{ neighbors}\\\\cell(i) \ne cell(j)}} J(\tau_i,\tau_j)
```

## Volume

The volume penalty simply sums the square-errors between a cell's current number of grid sites (i.e. volume) and its desired number of grids sites. 

```math
\text{Penalty} = (V - V_d)^2
```

In our case we set the desired volume was set to 20 in the `CellState` function. 

```@example CPMmodel
cpm.state.volumes[1:3]
```

And currently, all our cell volumes match our desired volume so no penalties are added. The value of 5 given to the `VolumePenalty` function is a constant that multiplies the total volume penalty. Larger values put more weight on this penalty. 

```math
\text{Total Volume Penalty} = 5 \sum_{i=1}^{3} (V_i - 20)^2
```

# Acceptance Probability

CPM uses a Metropolis–Hastings algorithm to determine if the update should be accepted. This is useful because unfavorable changes to the cells can still be accepted with a certain probability which will help avoid local minima. 

First we calculate the difference in energy $\Delta H$ between the original and proposed configurations. If the total energy decreases we automatically accept, otherwise we accept with a probability according to a Boltzmann distribution.

```math
P = e^{-\Delta H / T}
```

where $P$ is the probability of acceptance and $T$ is a "temperature" parameter. Larger temperatures coincide with less favorable events being accepted. 

# History and Criticism 

CPMs are based off of [Ising models](https://en.wikipedia.org/wiki/Ising_model) from statistical mechanics. However, a lot of physical meaning behind the parameters and ideas borrowed from these models has been lost. For example, the adhesion matrix $J$ used to represent atomic spin-spin interactions between magnetically charged atoms. The temperature parameter $T$ was literal temperature. This loss of physical meaning is one major criticism for CPMs which researchers are actively trying to remedy. 