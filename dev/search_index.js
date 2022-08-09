var documenterSearchIndex = {"docs":
[{"location":"ExampleGallery/LetsGetMoving/LetsGetMoving/","page":"Let's Get Moving","title":"Let's Get Moving","text":"EditURL = \"<unknown>/docs/src/ExampleGallery/LetsGetMoving/LetsGetMoving.jl\"","category":"page"},{"location":"ExampleGallery/LetsGetMoving/LetsGetMoving/#Let's-Get-Moving","page":"Let's Get Moving","title":"Let's Get Moving","text":"","category":"section"},{"location":"ExampleGallery/LetsGetMoving/LetsGetMoving/","page":"Let's Get Moving","title":"Let's Get Moving","text":"Many cells have the ability to move within their environment through the contraction of actin filaments. This mechanism leads to cells performing an \"intermittent random walk\" which is characterized by periods of persistent movement, followed by periods of being stationary.","category":"page"},{"location":"ExampleGallery/LetsGetMoving/LetsGetMoving/","page":"Let's Get Moving","title":"Let's Get Moving","text":"Let's see how we can add this kind of behavior to our model,","category":"page"},{"location":"ExampleGallery/LetsGetMoving/LetsGetMoving/","page":"Let's Get Moving","title":"Let's Get Moving","text":"Start by loading in the CellularPotts.jl package and creating a space where cells can exist.","category":"page"},{"location":"ExampleGallery/LetsGetMoving/LetsGetMoving/","page":"Let's Get Moving","title":"Let's Get Moving","text":"using CellularPotts\nspace = CellSpace(200,200)","category":"page"},{"location":"ExampleGallery/LetsGetMoving/LetsGetMoving/","page":"Let's Get Moving","title":"Let's Get Moving","text":"200×200 Periodic 8-Neighbor CellSpace{2,Int64}","category":"page"},{"location":"ExampleGallery/LetsGetMoving/LetsGetMoving/","page":"Let's Get Moving","title":"Let's Get Moving","text":"Much like in the HelloWorld.jl example, we create a single cell that averages 500 pixels in size.","category":"page"},{"location":"ExampleGallery/LetsGetMoving/LetsGetMoving/","page":"Let's Get Moving","title":"Let's Get Moving","text":"initialCellState = CellTable(\n    [:Epithelial],\n    [500],\n    [1]);","category":"page"},{"location":"ExampleGallery/LetsGetMoving/LetsGetMoving/","page":"Let's Get Moving","title":"Let's Get Moving","text":"The cell will be positioned at the halfway point within the space.","category":"page"},{"location":"ExampleGallery/LetsGetMoving/LetsGetMoving/","page":"Let's Get Moving","title":"Let's Get Moving","text":"positions = [size(space) .÷ 2]","category":"page"},{"location":"ExampleGallery/LetsGetMoving/LetsGetMoving/","page":"Let's Get Moving","title":"Let's Get Moving","text":"1-element Vector{Tuple{Int64, Int64}}:\n (100, 100)","category":"page"},{"location":"ExampleGallery/LetsGetMoving/LetsGetMoving/","page":"Let's Get Moving","title":"Let's Get Moving","text":"And that property is added to the CellTable","category":"page"},{"location":"ExampleGallery/LetsGetMoving/LetsGetMoving/","page":"Let's Get Moving","title":"Let's Get Moving","text":"initialCellState = addcellproperty(initialCellState, :positions, positions)","category":"page"},{"location":"ExampleGallery/LetsGetMoving/LetsGetMoving/","page":"Let's Get Moving","title":"Let's Get Moving","text":"┌────────────┬─────────┬─────────┬─────────┬────────────────┬────────────┬───────────────────┬─────────────────────┐\n│      names │ cellIDs │ typeIDs │ volumes │ desiredVolumes │ perimeters │ desiredPerimeters │           positions │\n│     Symbol │   Int64 │   Int64 │   Int64 │          Int64 │      Int64 │             Int64 │ Tuple{Int64, Int64} │\n├────────────┼─────────┼─────────┼─────────┼────────────────┼────────────┼───────────────────┼─────────────────────┤\n│     Medium │       0 │       0 │       0 │              0 │          0 │                 0 │          (100, 100) │\n│ Epithelial │       1 │       1 │       0 │            500 │          0 │               264 │          (100, 100) │\n└────────────┴─────────┴─────────┴─────────┴────────────────┴────────────┴───────────────────┴─────────────────────┘\n","category":"page"},{"location":"ExampleGallery/LetsGetMoving/LetsGetMoving/","page":"Let's Get Moving","title":"Let's Get Moving","text":"Now the important part. To enable this type of cellular movement, we can add a MigrationPenalty to the model. This penalty requires 3 inputs:","category":"page"},{"location":"ExampleGallery/LetsGetMoving/LetsGetMoving/","page":"Let's Get Moving","title":"Let's Get Moving","text":"A maximum activity\nA scaling factor for this penalty (one for each cell type)\nThe size of the space","category":"page"},{"location":"ExampleGallery/LetsGetMoving/LetsGetMoving/","page":"Let's Get Moving","title":"Let's Get Moving","text":"penalties = [\n    AdhesionPenalty([0 30;\n                    30 30]),\n    VolumePenalty([5]),\n    PerimeterPenalty([0,10]),\n    MigrationPenalty(50, [50], size(space))\n    ]","category":"page"},{"location":"ExampleGallery/LetsGetMoving/LetsGetMoving/","page":"Let's Get Moving","title":"Let's Get Moving","text":"4-element Vector{Penalty}:\n AdhesionPenalty([0 30; 30 30])\n VolumePenalty([0, 5])\n PerimeterPenalty([0, 0, 10], 0, 0)\n MigrationPenalty(50, [0, 50], sparse(Int64[], Int64[], Int64[], 200, 200))","category":"page"},{"location":"ExampleGallery/LetsGetMoving/LetsGetMoving/","page":"Let's Get Moving","title":"Let's Get Moving","text":"Create a Cellular Potts Model object","category":"page"},{"location":"ExampleGallery/LetsGetMoving/LetsGetMoving/","page":"Let's Get Moving","title":"Let's Get Moving","text":"cpm = CellPotts(space, initialCellState, penalties);","category":"page"},{"location":"ExampleGallery/LetsGetMoving/LetsGetMoving/","page":"Let's Get Moving","title":"Let's Get Moving","text":"We can adjust the \"temperature\" of the model (which defaults to 20). Higher temperatures will make model updates that increase the overall energy more likely. This is not necessary for this model, but a helpful feature to know.","category":"page"},{"location":"ExampleGallery/LetsGetMoving/LetsGetMoving/","page":"Let's Get Moving","title":"Let's Get Moving","text":"cpm.temperature = 25.0","category":"page"},{"location":"ExampleGallery/LetsGetMoving/LetsGetMoving/","page":"Let's Get Moving","title":"Let's Get Moving","text":"25.0","category":"page"},{"location":"ExampleGallery/LetsGetMoving/LetsGetMoving/","page":"Let's Get Moving","title":"Let's Get Moving","text":"Our cell still needs to be placed into the space. This can be done using the positionCellsRandom!() function or because we have a \"positions\" property, we can use the positionCells!() function.","category":"page"},{"location":"ExampleGallery/LetsGetMoving/LetsGetMoving/","page":"Let's Get Moving","title":"Let's Get Moving","text":"positionCells!(cpm)","category":"page"},{"location":"ExampleGallery/LetsGetMoving/LetsGetMoving/","page":"Let's Get Moving","title":"Let's Get Moving","text":"Our model is more ready for simulation! This can be done using the using the ModelStep! function, interactively through the CellGUI function, or recorded as a gif using recordCPM. Any options to the GLMakie record function can be passed through.","category":"page"},{"location":"ExampleGallery/LetsGetMoving/LetsGetMoving/","page":"Let's Get Moving","title":"Let's Get Moving","text":"recordCPM(\"LetsGetMoving.gif\", cpm)","category":"page"},{"location":"ExampleGallery/LetsGetMoving/LetsGetMoving/","page":"Let's Get Moving","title":"Let's Get Moving","text":"<p style=\"text-align:center;\">\n    <img title=\"\" src=\"https://github.com/RobertGregg/CellularPotts.jl/blob/master/docs/src/ExampleGallery/LetsGetMoving/LetsGetMoving.gif?raw=true\" width=\"445\">\n</p>","category":"page"},{"location":"ExampleGallery/LetsGetMoving/LetsGetMoving/","page":"Let's Get Moving","title":"Let's Get Moving","text":"","category":"page"},{"location":"ExampleGallery/LetsGetMoving/LetsGetMoving/","page":"Let's Get Moving","title":"Let's Get Moving","text":"This page was generated using Literate.jl.","category":"page"},{"location":"ExampleGallery/OnPatrol/OnPatrol/","page":"On Patrol","title":"On Patrol","text":"EditURL = \"<unknown>/docs/src/ExampleGallery/OnPatrol/OnPatrol.jl\"","category":"page"},{"location":"ExampleGallery/OnPatrol/OnPatrol/#On-Patrol","page":"On Patrol","title":"On Patrol","text":"","category":"section"},{"location":"ExampleGallery/OnPatrol/OnPatrol/","page":"On Patrol","title":"On Patrol","text":"Here we combine ideas from the HelloWorld.jl and LetsGetMoving.jl examples to simulate T-Cells patrolling through a dense layer of epithelial cells.","category":"page"},{"location":"ExampleGallery/OnPatrol/OnPatrol/","page":"On Patrol","title":"On Patrol","text":"Start by loading in the CellularPotts.jl package and creating a space where cells can exist.","category":"page"},{"location":"ExampleGallery/OnPatrol/OnPatrol/","page":"On Patrol","title":"On Patrol","text":"using CellularPotts\nspace = CellSpace(200,200)","category":"page"},{"location":"ExampleGallery/OnPatrol/OnPatrol/","page":"On Patrol","title":"On Patrol","text":"200×200 Periodic 8-Neighbor CellSpace{2,Int64}","category":"page"},{"location":"ExampleGallery/OnPatrol/OnPatrol/","page":"On Patrol","title":"On Patrol","text":"Initialize a new CellTable with 75 epithelial cells and 5 T-Cells","category":"page"},{"location":"ExampleGallery/OnPatrol/OnPatrol/","page":"On Patrol","title":"On Patrol","text":"initialCellState = CellTable(\n    [:Epithelial, :TCell],\n    [500, 400],\n    [75, 5])","category":"page"},{"location":"ExampleGallery/OnPatrol/OnPatrol/","page":"On Patrol","title":"On Patrol","text":"┌────────────┬─────────┬─────────┬─────────┬────────────────┬────────────┬───────────────────┐\n│      names │ cellIDs │ typeIDs │ volumes │ desiredVolumes │ perimeters │ desiredPerimeters │\n│     Symbol │   Int64 │   Int64 │   Int64 │          Int64 │      Int64 │             Int64 │\n├────────────┼─────────┼─────────┼─────────┼────────────────┼────────────┼───────────────────┤\n│     Medium │       0 │       0 │       0 │              0 │          0 │                 0 │\n│ Epithelial │       1 │       1 │       0 │            500 │          0 │               264 │\n│ Epithelial │       2 │       1 │       0 │            500 │          0 │               264 │\n│ Epithelial │       3 │       1 │       0 │            500 │          0 │               264 │\n│ Epithelial │       4 │       1 │       0 │            500 │          0 │               264 │\n│ Epithelial │       5 │       1 │       0 │            500 │          0 │               264 │\n│ Epithelial │       6 │       1 │       0 │            500 │          0 │               264 │\n│ Epithelial │       7 │       1 │       0 │            500 │          0 │               264 │\n│ Epithelial │       8 │       1 │       0 │            500 │          0 │               264 │\n│ Epithelial │       9 │       1 │       0 │            500 │          0 │               264 │\n│ Epithelial │      10 │       1 │       0 │            500 │          0 │               264 │\n│ Epithelial │      11 │       1 │       0 │            500 │          0 │               264 │\n│ Epithelial │      12 │       1 │       0 │            500 │          0 │               264 │\n│ Epithelial │      13 │       1 │       0 │            500 │          0 │               264 │\n│ Epithelial │      14 │       1 │       0 │            500 │          0 │               264 │\n│ Epithelial │      15 │       1 │       0 │            500 │          0 │               264 │\n│ Epithelial │      16 │       1 │       0 │            500 │          0 │               264 │\n│ Epithelial │      17 │       1 │       0 │            500 │          0 │               264 │\n│ Epithelial │      18 │       1 │       0 │            500 │          0 │               264 │\n│ Epithelial │      19 │       1 │       0 │            500 │          0 │               264 │\n│ Epithelial │      20 │       1 │       0 │            500 │          0 │               264 │\n│ Epithelial │      21 │       1 │       0 │            500 │          0 │               264 │\n│ Epithelial │      22 │       1 │       0 │            500 │          0 │               264 │\n│ Epithelial │      23 │       1 │       0 │            500 │          0 │               264 │\n│ Epithelial │      24 │       1 │       0 │            500 │          0 │               264 │\n│ Epithelial │      25 │       1 │       0 │            500 │          0 │               264 │\n│ Epithelial │      26 │       1 │       0 │            500 │          0 │               264 │\n│ Epithelial │      27 │       1 │       0 │            500 │          0 │               264 │\n│ Epithelial │      28 │       1 │       0 │            500 │          0 │               264 │\n│ Epithelial │      29 │       1 │       0 │            500 │          0 │               264 │\n│ Epithelial │      30 │       1 │       0 │            500 │          0 │               264 │\n│ Epithelial │      31 │       1 │       0 │            500 │          0 │               264 │\n│ Epithelial │      32 │       1 │       0 │            500 │          0 │               264 │\n│ Epithelial │      33 │       1 │       0 │            500 │          0 │               264 │\n│ Epithelial │      34 │       1 │       0 │            500 │          0 │               264 │\n│ Epithelial │      35 │       1 │       0 │            500 │          0 │               264 │\n│ Epithelial │      36 │       1 │       0 │            500 │          0 │               264 │\n│ Epithelial │      37 │       1 │       0 │            500 │          0 │               264 │\n│ Epithelial │      38 │       1 │       0 │            500 │          0 │               264 │\n│ Epithelial │      39 │       1 │       0 │            500 │          0 │               264 │\n│ Epithelial │      40 │       1 │       0 │            500 │          0 │               264 │\n│ Epithelial │      41 │       1 │       0 │            500 │          0 │               264 │\n│ Epithelial │      42 │       1 │       0 │            500 │          0 │               264 │\n│ Epithelial │      43 │       1 │       0 │            500 │          0 │               264 │\n│ Epithelial │      44 │       1 │       0 │            500 │          0 │               264 │\n│ Epithelial │      45 │       1 │       0 │            500 │          0 │               264 │\n│ Epithelial │      46 │       1 │       0 │            500 │          0 │               264 │\n│ Epithelial │      47 │       1 │       0 │            500 │          0 │               264 │\n│ Epithelial │      48 │       1 │       0 │            500 │          0 │               264 │\n│ Epithelial │      49 │       1 │       0 │            500 │          0 │               264 │\n│ Epithelial │      50 │       1 │       0 │            500 │          0 │               264 │\n│ Epithelial │      51 │       1 │       0 │            500 │          0 │               264 │\n│ Epithelial │      52 │       1 │       0 │            500 │          0 │               264 │\n│ Epithelial │      53 │       1 │       0 │            500 │          0 │               264 │\n│ Epithelial │      54 │       1 │       0 │            500 │          0 │               264 │\n│ Epithelial │      55 │       1 │       0 │            500 │          0 │               264 │\n│ Epithelial │      56 │       1 │       0 │            500 │          0 │               264 │\n│ Epithelial │      57 │       1 │       0 │            500 │          0 │               264 │\n│ Epithelial │      58 │       1 │       0 │            500 │          0 │               264 │\n│ Epithelial │      59 │       1 │       0 │            500 │          0 │               264 │\n│ Epithelial │      60 │       1 │       0 │            500 │          0 │               264 │\n│ Epithelial │      61 │       1 │       0 │            500 │          0 │               264 │\n│ Epithelial │      62 │       1 │       0 │            500 │          0 │               264 │\n│ Epithelial │      63 │       1 │       0 │            500 │          0 │               264 │\n│ Epithelial │      64 │       1 │       0 │            500 │          0 │               264 │\n│ Epithelial │      65 │       1 │       0 │            500 │          0 │               264 │\n│ Epithelial │      66 │       1 │       0 │            500 │          0 │               264 │\n│ Epithelial │      67 │       1 │       0 │            500 │          0 │               264 │\n│ Epithelial │      68 │       1 │       0 │            500 │          0 │               264 │\n│ Epithelial │      69 │       1 │       0 │            500 │          0 │               264 │\n│ Epithelial │      70 │       1 │       0 │            500 │          0 │               264 │\n│ Epithelial │      71 │       1 │       0 │            500 │          0 │               264 │\n│ Epithelial │      72 │       1 │       0 │            500 │          0 │               264 │\n│ Epithelial │      73 │       1 │       0 │            500 │          0 │               264 │\n│ Epithelial │      74 │       1 │       0 │            500 │          0 │               264 │\n│ Epithelial │      75 │       1 │       0 │            500 │          0 │               264 │\n│      TCell │      76 │       2 │       0 │            400 │          0 │               236 │\n│      TCell │      77 │       2 │       0 │            400 │          0 │               236 │\n│      TCell │      78 │       2 │       0 │            400 │          0 │               236 │\n│      TCell │      79 │       2 │       0 │            400 │          0 │               236 │\n│      TCell │      80 │       2 │       0 │            400 │          0 │               236 │\n└────────────┴─────────┴─────────┴─────────┴────────────────┴────────────┴───────────────────┘\n","category":"page"},{"location":"ExampleGallery/OnPatrol/OnPatrol/","page":"On Patrol","title":"On Patrol","text":"Note that for the MigrationPenalty we set the epithelial cell's scaling factor to zero. This effectively removes this penalty from the cell type.","category":"page"},{"location":"ExampleGallery/OnPatrol/OnPatrol/","page":"On Patrol","title":"On Patrol","text":"penalties = [\n    AdhesionPenalty([0 20 20;\n                    20 20 30;\n                    20 30 50]),\n    VolumePenalty([10,10]),\n    PerimeterPenalty([0,5]),\n    MigrationPenalty(50, [0,100], size(space))\n    ]","category":"page"},{"location":"ExampleGallery/OnPatrol/OnPatrol/","page":"On Patrol","title":"On Patrol","text":"4-element Vector{Penalty}:\n AdhesionPenalty([0 20 20; 20 20 30; 20 30 50])\n VolumePenalty([0, 10, 10])\n PerimeterPenalty([0, 0, 5], 0, 0)\n MigrationPenalty(50, [0, 0, 100], sparse(Int64[], Int64[], Int64[], 200, 200))","category":"page"},{"location":"ExampleGallery/OnPatrol/OnPatrol/","page":"On Patrol","title":"On Patrol","text":"Create a new CellPotts model.","category":"page"},{"location":"ExampleGallery/OnPatrol/OnPatrol/","page":"On Patrol","title":"On Patrol","text":"cpm = CellPotts(space, initialCellState, penalties)","category":"page"},{"location":"ExampleGallery/OnPatrol/OnPatrol/","page":"On Patrol","title":"On Patrol","text":"Cell Potts Model:\nGrid: 200×200\nCell Counts: [Epithelial → 75] [TCell → 5] [Total → 80]\nModel Penalties: Adhesion Migration Perimeter Volume\nTemperature: 20.0\nSteps: 0","category":"page"},{"location":"ExampleGallery/OnPatrol/OnPatrol/","page":"On Patrol","title":"On Patrol","text":"Here we did not specify the positions of the cells, so they be randomly added to the space.","category":"page"},{"location":"ExampleGallery/OnPatrol/OnPatrol/","page":"On Patrol","title":"On Patrol","text":"positionCellsRandom!(cpm)","category":"page"},{"location":"ExampleGallery/OnPatrol/OnPatrol/","page":"On Patrol","title":"On Patrol","text":"Record the simulation","category":"page"},{"location":"ExampleGallery/OnPatrol/OnPatrol/","page":"On Patrol","title":"On Patrol","text":"recordCPM(\"OnPatrol.gif\", cpm)","category":"page"},{"location":"ExampleGallery/OnPatrol/OnPatrol/","page":"On Patrol","title":"On Patrol","text":"<p style=\"text-align:center;\">\n    <img title=\"\" src=\"https://github.com/RobertGregg/CellularPotts.jl/blob/master/docs/src/ExampleGallery/OnPatrol/OnPatrol.gif?raw=true\" width=\"445\">\n</p>","category":"page"},{"location":"ExampleGallery/OnPatrol/OnPatrol/","page":"On Patrol","title":"On Patrol","text":"","category":"page"},{"location":"ExampleGallery/OnPatrol/OnPatrol/","page":"On Patrol","title":"On Patrol","text":"This page was generated using Literate.jl.","category":"page"},{"location":"ExampleGallery/HelloWorld/HelloWorld/","page":"Hello World","title":"Hello World","text":"EditURL = \"<unknown>/docs/src/ExampleGallery/HelloWorld/HelloWorld.jl\"","category":"page"},{"location":"ExampleGallery/HelloWorld/HelloWorld/#Hello-World","page":"Hello World","title":"Hello World","text":"","category":"section"},{"location":"ExampleGallery/HelloWorld/HelloWorld/","page":"Hello World","title":"Hello World","text":"In this example, we'll specify a single stationary cell in the center of the grid.","category":"page"},{"location":"ExampleGallery/HelloWorld/HelloWorld/","page":"Hello World","title":"Hello World","text":"We start by loading in the CellularPotts.jl package and creating a space where cells can exist.","category":"page"},{"location":"ExampleGallery/HelloWorld/HelloWorld/","page":"Hello World","title":"Hello World","text":"using CellularPotts\n\nspace = CellSpace(50,50; wrapAround=true, cellNeighbors=mooreNeighbors)","category":"page"},{"location":"ExampleGallery/HelloWorld/HelloWorld/","page":"Hello World","title":"Hello World","text":"50×50 Periodic 8-Neighbor CellSpace{2,Int64}","category":"page"},{"location":"ExampleGallery/HelloWorld/HelloWorld/","page":"Hello World","title":"Hello World","text":"Here we create a 50 by 50 square grid with periodic boundary conditions where grid locations are connected to their 8 closest neighbors (4-cell neighborhoods are also available using the vonNeumannNeighbors function). By default, wrapAround is set to true and cellNeighbors uses the 8 closest neighbors.","category":"page"},{"location":"ExampleGallery/HelloWorld/HelloWorld/","page":"Hello World","title":"Hello World","text":"Next we need to initialize a table of cell information to put into this space.","category":"page"},{"location":"ExampleGallery/HelloWorld/HelloWorld/","page":"Hello World","title":"Hello World","text":"initialCellState = CellTable(\n    [:Epithelial],\n    [500],\n    [1])","category":"page"},{"location":"ExampleGallery/HelloWorld/HelloWorld/","page":"Hello World","title":"Hello World","text":"┌────────────┬─────────┬─────────┬─────────┬────────────────┬────────────┬───────────────────┐\n│      names │ cellIDs │ typeIDs │ volumes │ desiredVolumes │ perimeters │ desiredPerimeters │\n│     Symbol │   Int64 │   Int64 │   Int64 │          Int64 │      Int64 │             Int64 │\n├────────────┼─────────┼─────────┼─────────┼────────────────┼────────────┼───────────────────┤\n│     Medium │       0 │       0 │       0 │              0 │          0 │                 0 │\n│ Epithelial │       1 │       1 │       0 │            500 │          0 │               264 │\n└────────────┴─────────┴─────────┴─────────┴────────────────┴────────────┴───────────────────┘\n","category":"page"},{"location":"ExampleGallery/HelloWorld/HelloWorld/","page":"Hello World","title":"Hello World","text":"The CellTable() function populates a table detailing the current cell state. The 3 required inputs are:","category":"page"},{"location":"ExampleGallery/HelloWorld/HelloWorld/","page":"Hello World","title":"Hello World","text":"A list of cell types\nA list of desired cell sizes (volumes)\nA list of cell counts for each cell type","category":"page"},{"location":"ExampleGallery/HelloWorld/HelloWorld/","page":"Hello World","title":"Hello World","text":"The inputs are simple in this case. We want one cell type called \"Epithelial\" with a size of 500 pixels and we want only one of them.","category":"page"},{"location":"ExampleGallery/HelloWorld/HelloWorld/","page":"Hello World","title":"Hello World","text":"The table CellTable() generates has each row representing a cell and each column listing a property given to that cell. Other information, like the column's type, is also provided.","category":"page"},{"location":"ExampleGallery/HelloWorld/HelloWorld/","page":"Hello World","title":"Hello World","text":"The first row will always show properties for \"Medium\", the name given to grid locations without a cell type. Most values related to Medium are either default or missing altogether. Here we see our one epithelial cell has a desired volume of 500 and perimeter of 264 which is the minimal perimeter penalty calculated from the desired volume.","category":"page"},{"location":"ExampleGallery/HelloWorld/HelloWorld/","page":"Hello World","title":"Hello World","text":"Additional properties can be added to our cells using the addcellproperty function. In this model we can provide a special property called positions with a single default value","category":"page"},{"location":"ExampleGallery/HelloWorld/HelloWorld/","page":"Hello World","title":"Hello World","text":"positions = [(25,25)]\n\ninitialCellState = addcellproperty(initialCellState, :positions, positions)","category":"page"},{"location":"ExampleGallery/HelloWorld/HelloWorld/","page":"Hello World","title":"Hello World","text":"┌────────────┬─────────┬─────────┬─────────┬────────────────┬────────────┬───────────────────┬─────────────────────┐\n│      names │ cellIDs │ typeIDs │ volumes │ desiredVolumes │ perimeters │ desiredPerimeters │           positions │\n│     Symbol │   Int64 │   Int64 │   Int64 │          Int64 │      Int64 │             Int64 │ Tuple{Int64, Int64} │\n├────────────┼─────────┼─────────┼─────────┼────────────────┼────────────┼───────────────────┼─────────────────────┤\n│     Medium │       0 │       0 │       0 │              0 │          0 │                 0 │            (25, 25) │\n│ Epithelial │       1 │       1 │       0 │            500 │          0 │               264 │            (25, 25) │\n└────────────┴─────────┴─────────┴─────────┴────────────────┴────────────┴───────────────────┴─────────────────────┘\n","category":"page"},{"location":"ExampleGallery/HelloWorld/HelloWorld/","page":"Hello World","title":"Hello World","text":"Looking at our updated table, we can see the newly added property.","category":"page"},{"location":"ExampleGallery/HelloWorld/HelloWorld/","page":"Hello World","title":"Hello World","text":"Now that we have a space and a cell to fill it with, we need to provide a list of model penalties. A number of default penalties exist and you can even create your own custom penalties. Here we only include an AdhesionPenalty which encourages grid locations with the same cell type to stick together and a VolumePenalty which penalizes cells that deviate from their desired volume.","category":"page"},{"location":"ExampleGallery/HelloWorld/HelloWorld/","page":"Hello World","title":"Hello World","text":"penalties = [\n    AdhesionPenalty([0 20;\n                     20 0]),\n    VolumePenalty([5])\n    ]","category":"page"},{"location":"ExampleGallery/HelloWorld/HelloWorld/","page":"Hello World","title":"Hello World","text":"2-element Vector{Penalty}:\n AdhesionPenalty([0 20; 20 0])\n VolumePenalty([0, 5])","category":"page"},{"location":"ExampleGallery/HelloWorld/HelloWorld/","page":"Hello World","title":"Hello World","text":"AdhesionPenalty requires a symmetric matrix J where J[n,m] gives the adhesion penalty for cells with types n and m. In this model we penalize Epithelial cell locations adjacent to Medium. The VolumePenalty needs a vector of scaling factors (one for each cell type) that either increase or decrease the volume penalty contribution to the overall penalty. The scaling factor for :Medium is automatically set to zero.","category":"page"},{"location":"ExampleGallery/HelloWorld/HelloWorld/","page":"Hello World","title":"Hello World","text":"Now we can take these three objects and create a Cellular Potts Model object.","category":"page"},{"location":"ExampleGallery/HelloWorld/HelloWorld/","page":"Hello World","title":"Hello World","text":"cpm = CellPotts(space, initialCellState, penalties)","category":"page"},{"location":"ExampleGallery/HelloWorld/HelloWorld/","page":"Hello World","title":"Hello World","text":"Cell Potts Model:\nGrid: 50×50\nCell Counts: [Epithelial → 1] [Total → 1]\nModel Penalties: Adhesion Volume\nTemperature: 20.0\nSteps: 0","category":"page"},{"location":"ExampleGallery/HelloWorld/HelloWorld/","page":"Hello World","title":"Hello World","text":"Calling this object gives a quick summary of the model's current state. Note that a \"temperature\" of 20 is given to the model by default. Higher temperatures allow the model to more likely accept changes that increase the overall penalty (e.g. cells could deviate further from their desired volume). The model object also tracks how many time steps have been performed.","category":"page"},{"location":"ExampleGallery/HelloWorld/HelloWorld/","page":"Hello World","title":"Hello World","text":"Our cell still needs to be placed into the space. This can be done using the positionCellsRandom!() function or because we have a \"positions\" property, we can use the positionCells!() function.","category":"page"},{"location":"ExampleGallery/HelloWorld/HelloWorld/","page":"Hello World","title":"Hello World","text":"positionCells!(cpm)","category":"page"},{"location":"ExampleGallery/HelloWorld/HelloWorld/","page":"Hello World","title":"Hello World","text":"Our model is more ready for simulation! This can be done using the using the ModelStep! function, interactively through the CellGUI function, or recorded as a gif using recordCPM","category":"page"},{"location":"ExampleGallery/HelloWorld/HelloWorld/","page":"Hello World","title":"Hello World","text":"recordCPM(\"HelloWorld.gif\", cpm)","category":"page"},{"location":"ExampleGallery/HelloWorld/HelloWorld/","page":"Hello World","title":"Hello World","text":"<p style=\"text-align:center;\">\n    <img title=\"\" src=\"https://github.com/RobertGregg/CellularPotts.jl/blob/master/docs/src/ExampleGallery/HelloWorld/HelloWorld.gif?raw=true\" width=\"445\">\n</p>","category":"page"},{"location":"ExampleGallery/HelloWorld/HelloWorld/","page":"Hello World","title":"Hello World","text":"","category":"page"},{"location":"ExampleGallery/HelloWorld/HelloWorld/","page":"Hello World","title":"Hello World","text":"This page was generated using Literate.jl.","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"DocTestSetup = quote\n    using CellularPotts\nend","category":"page"},{"location":"#Introduction","page":"Introduction","title":"Introduction","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"CPMs simulate a collection of cells interacting with one another. These interactions can range anywhere from simple contact between cells to complex cytokine communication.","category":"page"},{"location":"#Installation","page":"Introduction","title":"Installation","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"using Pkg\nPkg.add(\"CellularPotts\")","category":"page"},{"location":"#Example-Gallery","page":"Introduction","title":"Example Gallery","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"tip: Tip\nClick on an example to go to example page","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"\n<style>\n    .Gallery {\n        display: grid;\n        grid-template-columns: repeat(3,auto);\n        grid-gap: 1%;\n    }  \n\n    h3 {\n        text-align: center;\n    }\n</style>\n\n<div class=\"Gallery\">\n\n    <div>\n        <h3>Hello world</h3>\n        <a href=\"https://robertgregg.github.io/CellularPotts.jl/dev/ExampleGallery/HelloWorld/HelloWorld/\">\n            <img title=\"\" src=\"https://github.com/RobertGregg/CellularPotts.jl/blob/master/docs/src/ExampleGallery/HelloWorld/HelloWorld.gif?raw=true\">\n        </a>\n    </div>\n\n    <div>\n        <h3>Let's get moving</h3>\n        <a href=\"https://robertgregg.github.io/CellularPotts.jl/dev/ExampleGallery/LetsGetMoving/LetsGetMoving/\">\n            <img title=\"\" src=\"https://github.com/RobertGregg/CellularPotts.jl/blob/master/docs/src/ExampleGallery/LetsGetMoving/LetsGetMoving.gif?raw=true\">\n        </a>\n    </div>\n\n    <div>\n        <h3>On patrol</h3>\n        <a href=\"https://robertgregg.github.io/CellularPotts.jl/dev/ExampleGallery/OnPatrol/OnPatrol/\">\n            <img title=\"\" src=\"https://github.com/RobertGregg/CellularPotts.jl/blob/master/docs/src/ExampleGallery/OnPatrol/OnPatrol.gif?raw=true\">\n        </a>\n    </div>\n\n</div>\n","category":"page"},{"location":"API/#API","page":"API","title":"API","text":"","category":"section"},{"location":"API/","page":"API","title":"API","text":"CellularPotts.jl includes the following core functions.","category":"page"},{"location":"API/#Index","page":"API","title":"Index","text":"","category":"section"},{"location":"API/","page":"API","title":"API","text":"Pages = [\"API.md\"]","category":"page"},{"location":"API/#Full-docs","page":"API","title":"Full docs","text":"","category":"section"},{"location":"API/","page":"API","title":"API","text":"Modules = [CellularPotts]\nPages = [\"Core.jl\"]","category":"page"},{"location":"API/#CellularPotts.AdhesionPenalty","page":"API","title":"CellularPotts.AdhesionPenalty","text":"AdhesionPenalty(J::Matrix{Int})\n\nAn concrete type that penalizes neighboring grid locations from different cells.\n\nRequires a symmetric matrix J where J[n,m] gives the adhesion penality for cells with types n and m. J is zero-indexed meaning J[0,1] and J[1,0] corresponds to the :Medium ↔ :Cell1 adhesion penalty.\n\nNote: J is automatically transformed to be a zero-indexed offset array.\n\n\n\n\n\n","category":"type"},{"location":"API/#CellularPotts.CellPotts","page":"API","title":"CellularPotts.CellPotts","text":"CellPotts(space, initialCellState, penalties)\n\nA data container that holds information to run the cellular potts simulation.\n\nRequires three inputs:\n\nspace: a region where cells can exist, generated using CellSpace().\ninitialCellState: a table where rows are cells and columns are cell properties, generated using CellTable().\npenalties: a vector of penalties to append to the model.\n\n\n\n\n\n","category":"type"},{"location":"API/#CellularPotts.MigrationPenalty","page":"API","title":"CellularPotts.MigrationPenalty","text":"MigrationPenalty(maxAct, λ, gridSize)\n\nAn concrete type that encourages cells to protude and drag themselves forward.\n\nTwo integer parameters control how cells protude:\n\nmaxAct: A maximum activity a grid location can have\nλ: A parameter that controls the strength of this penalty\n\nIncreasing maxAct will cause grid locations to more likely protrude. Increasing λ will cause those protusions to reach farther away. \n\nMigrationPenalty also requires a list of cell types to apply the penalty to and the grid size (Space.gridSize).\n\n\n\n\n\n","category":"type"},{"location":"API/#CellularPotts.Penalty","page":"API","title":"CellularPotts.Penalty","text":"Penalty\n\nAn abstract type representing a constraint imposed onto the cellular potts model.\n\nTo add a new penalty, a new struct subtyping Penalty needs to be defined and the addPenalty!() function needs to be extended to include the new penalty.\n\nNote: variables associated with a new penalty may need to be offset such that index 0 maps to :Medium, index 1 maps to :Cell1, etc.\n\n\n\n\n\n","category":"type"},{"location":"API/#CellularPotts.PerimeterPenalty","page":"API","title":"CellularPotts.PerimeterPenalty","text":"PerimeterPenalty(λᵥ::Vector{Int})\n\nAn concrete type that penalizes cells that deviate from their desired perimeter.\n\nRequires a vector λₚ with n penalties where n is the number of cell types. λₚ is zero-indexed meaning λₚ[0] corresponds to the :Medium perimeter penalty (which is set to zero).\n\nNote: λₚ is automatically transformed to be a zero-indexed offset array and does not require the perimeter penalty for :Medium.\n\n\n\n\n\n","category":"type"},{"location":"API/#CellularPotts.VolumePenalty","page":"API","title":"CellularPotts.VolumePenalty","text":"VolumePenalty(λᵥ::Vector{Int})\n\nAn concrete type that penalizes cells that deviate from their desired volume.\n\nRequires a vector λᵥ with n penalties where n is the number of cell types. λᵥ is zero-indexed meaning λᵥ[0] corresponds to the :Medium volume penalty (which is set to zero).\n\nNote: λᵥ is automatically transformed to be a zero-indexed offset array and does not require the volume penalty for :Medium.\n\n\n\n\n\n","category":"type"},{"location":"API/#CellularPotts.countcells-Tuple{CellPotts}","page":"API","title":"CellularPotts.countcells","text":"countcells(cpm::CellPotts)\ncountcells(df::CellTable)\n\nCount the number of cells in the model \n\n\n\n\n\n","category":"method"},{"location":"API/#CellularPotts.countcelltypes-Tuple{CellPotts}","page":"API","title":"CellularPotts.countcelltypes","text":"countcelltypes(cpm::CellPotts)\ncountcelltypes(df::CellTable)\n\nCount the number of cell types in the model \n\n\n\n\n\n","category":"method"}]
}
