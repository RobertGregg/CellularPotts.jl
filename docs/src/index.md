```@meta
DocTestSetup = quote
    using CellularPotts
end
```

# Introduction

CPMs simulate a collection of cells interacting with one another. These interactions can range anywhere from simple contact between cells to complex cytokine communication.

# Installation

```julia
using Pkg
Pkg.add("CellularPotts")
```

# Simple example

```@raw html
<style>
    .Gallery {
        display: grid;
        grid-template-columns: repeat(4,auto);
        grid-gap: 1%;
    }  
</style>

<div class="Gallery">

    <div>
        <a href="https://robertgregg.github.io/CellularPotts.jl/dev/ExampleGallery/HelloWorld/HelloWorld/">
            <img title="" src="https://github.com/RobertGregg/CellularPotts.jl/blob/master/docs/src/ExampleGallery/HelloWorld/HelloWorld.gif?raw=true">
        </a>
    </div>

    <div>
        <a href="https://robertgregg.github.io/CellularPotts.jl/dev/ExampleGallery/HelloWorld/HelloWorld/">
            <img title="" src="https://github.com/RobertGregg/CellularPotts.jl/blob/master/docs/src/ExampleGallery/HelloWorld/HelloWorld.gif?raw=true">
        </a>
    </div>

    <div>
        <a href="https://robertgregg.github.io/CellularPotts.jl/dev/ExampleGallery/HelloWorld/HelloWorld/">
            <img title="" src="https://github.com/RobertGregg/CellularPotts.jl/blob/master/docs/src/ExampleGallery/HelloWorld/HelloWorld.gif?raw=true">
        </a>
    </div>

    <div>
        <a href="https://robertgregg.github.io/CellularPotts.jl/dev/ExampleGallery/HelloWorld/HelloWorld/">
            <img title="" src="https://github.com/RobertGregg/CellularPotts.jl/blob/master/docs/src/ExampleGallery/HelloWorld/HelloWorld.gif?raw=true">
        </a>
    </div>

    <div>
        <a href="https://robertgregg.github.io/CellularPotts.jl/dev/ExampleGallery/HelloWorld/HelloWorld/">
            <img title="" src="https://github.com/RobertGregg/CellularPotts.jl/blob/master/docs/src/ExampleGallery/HelloWorld/HelloWorld.gif?raw=true">
        </a>
    </div>

</div>

```

In this example, we'll specify a single stationary cell in the center of the grid. 

We start by loading in the `CellularPotts.jl` package and creating a space where cells can exist.

```jldoctest simpleExample
julia> using CellularPotts

julia> space = CellSpace(50,50; wrapAround=true, cellNeighbors=mooreNeighbors)
{2500, 10000} undirected simple Int64 graph
```

Here we create a 50 by 50 square grid with periodic boundary conditions where grid locations are connected to their 8 closest neighbors (4-cell neighborhoods are also available using the `vonNeumannNeighbors` function). By default, `wrapAround` is set to true and `cellNeighbors` uses the 8 closest neighbors. 

Next we need to initialize a table of cell information to put into this space.

```jldoctest simpleExample
julia> initialCellState = newCellState([:Epithelial], [500], [1])
┌────────────┬─────────┬─────────┬─────────┬────────────────┬────────────┬───────────────────┐
│      names │ cellIDs │ typeIDs │ volumes │ desiredVolumes │ perimeters │ desiredPerimeters │
│     Symbol │   Int64 │   Int64 │   Int64 │          Int64 │      Int64 │             Int64 │
├────────────┼─────────┼─────────┼─────────┼────────────────┼────────────┼───────────────────┤
│     Medium │       0 │       0 │       0 │              0 │          0 │                 0 │
│ Epithelial │       1 │       1 │       0 │            500 │          0 │               264 │
└────────────┴─────────┴─────────┴─────────┴────────────────┴────────────┴───────────────────┘
```

The `newCellState()` function populates a table detailing the current cell state. The 3 required inputs are:

1. A list of cell types

2. A list of desired cell sizes (volumes)

3. A list of cell counts for each cell type

The inputs are simple in this case. We want one cell type called "Epithelial" with a size of 500 pixels and we want only one of them.

The table `newCellState()` generates has each row representing a cell and each column listing a property given to that cell. Other information, like the column's type, is also provided.

The first row will always show properties for "Medium", the name given to grid locations without a cell type. Most values related to Medium are either default or missing altogether. Here we see our one epithelial cell has a desired volume of 500 and perimeter of 264 which is the minimal perimeter penalty calculated from the desired volume. 

Additional properties can be added to our cells. In this model we can provide a property called positions with a single default value

```jldoctest simpleExample
julia> positions = [(25,25)];

julia> initialCellState = addcellproperty(initialCellState, :positions, positions)
┌────────────┬─────────┬─────────┬─────────┬────────────────┬────────────┬───────────────────┬─────────────────────┐
│      names │ cellIDs │ typeIDs │ volumes │ desiredVolumes │ perimeters │ desiredPerimeters │           positions │
│     Symbol │   Int64 │   Int64 │   Int64 │          Int64 │      Int64 │             Int64 │ Tuple{Int64, Int64} │
├────────────┼─────────┼─────────┼─────────┼────────────────┼────────────┼───────────────────┼─────────────────────┤
│     Medium │       0 │       0 │       0 │              0 │          0 │                 0 │            (25, 25) │
│ Epithelial │       1 │       1 │       0 │            500 │          0 │               264 │            (25, 25) │
└────────────┴─────────┴─────────┴─────────┴────────────────┴────────────┴───────────────────┴─────────────────────┘
```

Looking at our updated table, we can see the newly added property.

Now that we have a space and a cell to fill it with, we need to provide a list of model penalties. A number of default penalties exist and you can even create your own custom penalties. Here we only include an `AdhesionPenalty` which encourages grid locations with the same cell type to stick together and a `VolumePenalty` which penalizes cells that deviate from their desired volume.

```jldoctest simpleExample
julia> penalties = [
           AdhesionPenalty([0 20;
                           20 0]),
           VolumePenalty([5])
           ]
2-element Vector{Penalty}:
 AdhesionPenalty([0 20; 20 0])
 VolumePenalty([0, 5])
```

`AdhesionPenalty` requires a symmetric matrix `J` where `J[n,m]` gives the adhesion penalty for cells with types n and m. In this model we penalize Epithelial cell locations adjacent to Medium. The `VolumePenalty` needs a vector of scaling factors (one for each cell type) that either increase or decrease the volume penalty contribution to the overall penalty. 

Now we can take these three objects and create a Cellular Potts Model object.

```jldoctest simpleExample
julia> cpm = CellPotts(space, initialCellState, penalties)
Cell Potts Model:
Grid: 50×50
Cell Counts: [Epithelial → 1] [Total → 1]
Model Penalties: Adhesion Volume
Temperature: 20.0
Steps: 0 
```

Calling this object gives a quick summary of the model's current state. Note that a "temperature" of 20 is given to the model by default. Higher temperatures allow the model to more likely accept changes that increase the overall penalty (e.g. cells could deviate further from their desired volume). The model object also tracks how many time steps have been performed. 

Our cell still needs to be placed into the space. This can be done using the `positionCellsRandom!()` function or because we have a "positions" property, we can use the `positionCells!()` function.

```julia
positionCells!(cpm)
```

Our model is more ready for simulation! This can be done using the using the `ModelStep!` function, interactively through the `CellGUI` function, or recorded as a gif using `recordCPM`

```julia
recordCPM("HelloWorld.gif", cpm)
```