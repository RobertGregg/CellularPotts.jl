# # Let's Get Moving

# Many cells have the ability to move within their environment through the contraction of actin filaments. This mechanism leads to cells performing an **"intermittent random walk"** which is characterized by periods of persistent movement, followed by periods of being stationary. 

# Let's see how we can add this kind of behavior to our model,

# Start by loading in the `CellularPotts.jl` package and creating a space where cells can exist.
using CellularPotts
space = CellSpace(100,100)

# Much like in the [HelloWorld.jl](https://robertgregg.github.io/CellularPotts.jl/dev/ExampleGallery/HelloWorld/HelloWorld/) example, we create a single cell that averages 500 pixels in size.

initialCellState = CellState(:Epithelial, 500, 1);

# The cell will be positioned at the halfway point within the space. 

positions = [size(space) .÷ 2]

# And that property is added to the CellState

initialCellState = addcellproperty(initialCellState, :positions, positions)

# Now the important part. To enable this type of cellular movement, we can add a `MigrationPenalty` to the model. This penalty requires 3 inputs:

# 1. A maximum activity 
# 2. A scaling factor for this penalty (one for each cell type)
# 3. The size of the space

penalties = [
    AdhesionPenalty([0 30;
                    30 30]),
    VolumePenalty([5]),
    PerimeterPenalty([0,10]),
    MigrationPenalty(50, [50], size(space))
    ]


# Create a Cellular Potts Model object

cpm = CellPotts(space, initialCellState, penalties);

# We can adjust the "temperature" of the model (which defaults to 20). Higher temperatures will make model updates that increase the overall energy more likely. This is not necessary for this model, but a helpful feature to know. 

cpm.temperature = 25.0

# Our model is more ready for simulation! This can be done using the using the `ModelStep!` function, interactively through the `CellGUI` function, or recorded as a gif using `recordCPM`. Any options to the GLMakie `record` function can be passed through.

recordCPM("LetsGetMoving.gif", cpm)