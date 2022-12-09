# # Hello World

# In this example, we'll specify a single stationary cell in the center of the grid. 

# We start by loading in the `CellularPotts.jl` package and creating a space where cells can exist.

using CellularPotts

space = CellSpace(50,50; isPeriodic=true, neighborhood=:moore)

# Here we create a 50 by 50 square grid with periodic boundary conditions where grid locations are connected to their 8 closest neighbors (4-cell neighborhoods are also available using the `:vonNeumann` option). By default, `isPeriodic` is set to true and `neighborhood` uses the 8 closest neighbors.

# Next we need to initialize a table of cell information to put into this space.

initialCellState = CellTable(
    [:Epithelial],
    [500],
    [1])

# The `CellTable()` function populates a table detailing the current cell state. The 3 required inputs are:

# 1. A list of cell types    
# 2. A list of desired cell sizes (volumes)
# 3. A list of cell counts for each cell type
    
# The inputs are simple in this case. We want one cell type called "Epithelial" with a size of 500 pixels and we want only one of them.
    
# The table `CellTable()` generates has each row representing a cell and each column listing a property given to that cell. Other information, like the column's type, is also provided.
    
# The first row will always show properties for "Medium", the name given to grid locations without a cell type. Most values related to Medium are either default or missing altogether. Here we see our one epithelial cell has a desired volume of 500 and perimeter of 264 which is the minimal perimeter penalty calculated from the desired volume. 
    
# Additional properties can be added to our cells using the `addcellproperty` function. In this model we can provide a special property called positions with a single default value

positions = [(25,25)]

initialCellState = addcellproperty(initialCellState, :positions, positions)

# Looking at our updated table, we can see the newly added property.

# Now that we have a space and a cell to fill it with, we need to provide a list of model penalties. A number of default penalties exist and you can even create your own custom penalties. Here we only include an `AdhesionPenalty` which encourages grid locations with the same cell type to stick together and a `VolumePenalty` which penalizes cells that deviate from their desired volume.

penalties = [
    AdhesionPenalty([0 20;
                     20 0]),
    VolumePenalty([5])
    ]

# `AdhesionPenalty` requires a symmetric matrix `J` where `J[n,m]` gives the adhesion penalty for cells with types n and m. In this model we penalize Epithelial cell locations adjacent to Medium. The `VolumePenalty` needs a vector of scaling factors (one for each cell type) that either increase or decrease the volume penalty contribution to the overall penalty. The scaling factor for `:Medium` is automatically set to zero.

# Now we can take these three objects and create a Cellular Potts Model object.

cpm = CellPotts(space, initialCellState, penalties)

# Calling this object gives a quick summary of the model's current state. Note that a "temperature" of 20 is given to the model by default. Higher temperatures allow the model to more likely accept changes that increase the overall penalty (e.g. cells could deviate further from their desired volume). The model object also tracks how many time steps have been performed. 

# !!! note 
# Because we added the `:positions` property to our CellTable, the CellPotts initializer placed the cells in the grid centered at those positions. If we did not specify a :positions property, cells would be randomly placed in the space.

# Our model is more ready for simulation! This can be done using the using the `ModelStep!` function, interactively through the `CellGUI` function, or recorded as a gif using `recordCPM`

recordCPM("HelloWorld.gif", cpm)