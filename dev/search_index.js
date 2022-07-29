var documenterSearchIndex = {"docs":
[{"location":"#Introduction","page":"Introduction","title":"Introduction","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"CPMs simulate a collection of cells interacting with one another. These interactions can range anywhere from simple contact between cells to complex cytokine communication.","category":"page"},{"location":"#Installation","page":"Introduction","title":"Installation","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"pkg> add CellularPotts","category":"page"},{"location":"#Simple-example","page":"Introduction","title":"Simple example","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"<img title=\"\" src=\"https://github.com/RobertGregg/CellularPotts.jl/blob/master/src/Gallery/HelloWorld/HelloWorld.gif\" alt=\"\" width=\"445\">","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"(Image: )","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"<img title=\"\" src=\"https://github.com/RobertGregg/CellularPotts.jl/blob/master/src/Gallery/HelloWorld/HelloWorld.gif?raw=true\" alt=\"\" width=\"445\">","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"(Image: )","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"In this example, we'll specify a single stationary cell in the center of the grid. ","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"We start by loading in the CellularPotts.jl package and creating a space where cells can exist.","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"using CellularPotts\n\nspace = CellSpace(50,50; wrapAround=true, cellNeighbors=mooreNeighbors)","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"Here we create a 50 by 50 square grid with periodic boundary conditions where grid locations are connected to their 8 closest neighbors (4-cell neighborhoods are also available using the vonNeumannNeighbors function). By default, wrapAround is set to true and cellNeighbors uses the 8 closest neighbors. ","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"Next we need to initialize a table of cell information to put into this space.","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"initialCellState = newCellState([:Epithelial], [500], [1]); ","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"The newCellState() function populates a table detailing the current cell state. The 3 required inputs are:","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"A list of cell types\nA list of desired cell sizes (volumes)\nA list of cell counts for each cell type","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"The inputs are simple in this case. We want one cell type called \"Epithelial\" with a size of 500 pixels and we want only one of them.","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"The table newCellState() generates has each row representing a cell and each column listing a property given to that cell. Other information, like the column's type, is also provided.","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"julia> initialCellState\n┌────────────┬─────────┬─────────┬─────────┬────────────────┬────────────┬───────────────────┐\n│      names │ cellIDs │ typeIDs │ volumes │ desiredVolumes │ perimeters │ desiredPerimeters │\n│     Symbol │   Int64 │   Int64 │   Int64 │          Int64 │      Int64 │             Int64 │\n├────────────┼─────────┼─────────┼─────────┼────────────────┼────────────┼───────────────────┤\n│     Medium │       0 │       0 │       0 │              0 │          0 │                 0 │\n│ Epithelial │       1 │       1 │       0 │            500 │          0 │               264 │\n└────────────┴─────────┴─────────┴─────────┴────────────────┴────────────┴───────────────────┘","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"The first row will always show properties for \"Medium\", the name given to grid locations without a cell type. Most values related to Medium are either default or missing altogether. Here we see our one epithelial cell has a desired volume of 500 and perimeter of 264 which is the minimal perimeter penalty calculated from the desired volume. ","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"Additional properties can be added to our cells. In this model we can provide a property called positions with a single default value","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"positions = [(25,25)]\ninitialCellState = addCellProperty(initialCellState, :positions, positions)","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"Looking at our updated table, we can see the newly added property.","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"julia> initialCellState\n┌────────────┬─────────┬─────────┬─────────┬────────────────┬────────────┬───────────────────┬─────────────────────┐\n│      names │ cellIDs │ typeIDs │ volumes │ desiredVolumes │ perimeters │ desiredPerimeters │           positions │\n│     Symbol │   Int64 │   Int64 │   Int64 │          Int64 │      Int64 │             Int64 │ Tuple{Int64, Int64} │\n├────────────┼─────────┼─────────┼─────────┼────────────────┼────────────┼───────────────────┼─────────────────────┤\n│     Medium │       0 │       0 │       0 │              0 │          0 │                 0 │            (25, 25) │\n│ Epithelial │       1 │       1 │       0 │            500 │          0 │               264 │            (25, 25) │\n└────────────┴─────────┴─────────┴─────────┴────────────────┴────────────┴───────────────────┴─────────────────────┘","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"Now that we have a space and a cell to fill it with, we need to provide a list of model penalties. A number of default penalties exist and you can even create your own custom penalties. Here we only include an AdhesionPenalty which encourages grid locations with the same cell type to stick together and a VolumePenalty which penalizes cells that deviate from their desired volume.","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"penalties = [\n    AdhesionPenalty([ 0  20;\n                     20  0]),\n    VolumePenalty([5])\n    ]","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"AdhesionPenalty requires a symmetric matrix J where J[n,m] gives the adhesion penalty for cells with types n and m. In this model we penalize Epithelial cell locations adjacent to Medium. The VolumePenalty needs a vector of scaling factors (one for each cell type) that either increase or decrease the volume penalty contribution to the overall penalty. ","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"Now we can take these three objects and create a Cellular Potts Model object.","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"julia> cpm = CellPotts(space, initialCellState, penalties);\n\njulia> cpm\nCell Potts Model:\nGrid: 50×50\nCell Counts: [Epithelial → 1] [Total → 1]\nModel Penalties: Adhesion Volume\nTemperature: 20.0\nSteps: 0 ","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"Calling this object gives a quick summary of the model's current state. Note that a \"temperature\" of 20 is given to the model by default. Higher temperatures allow the model to more likely accept changes that increase the overall penalty (e.g. cells could deviate further from their desired volume). The model object also tracks how many time steps have been performed. ","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"Our cell still needs to be placed into the space. This can be done using the positionCellsRandom!() function or because we have a \"positions\" property, we can use the positionCells!() function.","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"positionCells!(cpm)","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"Our model is more ready for simulation! This can be done using the using the ModelStep! function, interactively through the CellGUI() function, or recorded as a gif using cpmSaveGif","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"cpmSaveGif()","category":"page"},{"location":"API/#API","page":"API","title":"API","text":"","category":"section"},{"location":"API/","page":"API","title":"API","text":"CellularPotts.jl includes the following core functions.","category":"page"},{"location":"API/#Index","page":"API","title":"Index","text":"","category":"section"},{"location":"API/","page":"API","title":"API","text":"Pages = [\"API.md\"]","category":"page"},{"location":"API/#Full-docs","page":"API","title":"Full docs","text":"","category":"section"},{"location":"API/","page":"API","title":"API","text":"Modules = [CellularPotts]\nPages = [\"Core.jl\"]","category":"page"},{"location":"API/#CellularPotts.AdhesionPenalty","page":"API","title":"CellularPotts.AdhesionPenalty","text":"AdhesionPenalty(J::Matrix{Int})\n\nAn concrete type that penalizes neighboring grid locations from different cells.\n\nRequires a symmetric matrix J where J[n,m] gives the adhesion penality for cells with types n and m. J is zero-indexed meaning J[0,1] and J[1,0] corresponds to the :Medium ↔ :Cell1 adhesion penalty.\n\nNote: J is automatically transformed to be a zero-indexed offset array.\n\n\n\n\n\n","category":"type"},{"location":"API/#CellularPotts.CellPotts","page":"API","title":"CellularPotts.CellPotts","text":"CellPotts(space, initialCellState, penalties)\n\nA data container that holds information to run the cellular potts simulation.\n\nRequires three inputs:\n\nspace: a region where cells can exist, generated using CellSpace().\ninitialCellState: a table where rows are cells and columns are cell properties, generated using newCellState().\npenalties: a vector of penalties to append to the model.\n\n\n\n\n\n","category":"type"},{"location":"API/#CellularPotts.CellTable","page":"API","title":"CellularPotts.CellTable","text":"CellTable\n\nAn concrete type that stores a table where each row is a cell and each column is a cell property.\n\nThis table can be genereated using the newCellState() function.\n\n\n\n\n\n","category":"type"},{"location":"API/#CellularPotts.MigrationPenalty","page":"API","title":"CellularPotts.MigrationPenalty","text":"MigrationPenalty(maxAct, λ, gridSize)\n\nAn concrete type that encourages cells to protude and drag themselves forward.\n\nTwo integer parameters control how cells protude:\n\nmaxAct: A maximum activity a grid location can have\nλ: A parameter that controls the strength of this penalty\n\nIncreasing maxAct will cause grid locations to more likely protrude. Increasing λ will cause those protusions to reach farther away. \n\nMigrationPenalty also requires a list of cell types to apply the penalty to and the grid size (Space.gridSize).\n\n\n\n\n\n","category":"type"},{"location":"API/#CellularPotts.Penalty","page":"API","title":"CellularPotts.Penalty","text":"Penalty\n\nAn abstract type representing a constraint imposed onto the cellular potts model.\n\nTo add a new penalty, a new struct subtyping Penalty needs to be defined and the addPenalty!() function needs to be extended to include the new penalty.\n\nNote: variables associated with a new penalty may need to be offset such that index 0 maps to :Medium, index 1 maps to :Cell1, etc.\n\n\n\n\n\n","category":"type"},{"location":"API/#CellularPotts.PerimeterPenalty","page":"API","title":"CellularPotts.PerimeterPenalty","text":"PerimeterPenalty(λᵥ::Vector{Int})\n\nAn concrete type that penalizes cells that deviate from their desired perimeter.\n\nRequires a vector λₚ with n penalties where n is the number of cell types. λₚ is zero-indexed meaning λₚ[0] corresponds to the :Medium perimeter penalty (which is set to zero).\n\nNote: λₚ is automatically transformed to be a zero-indexed offset array and does not require the perimeter penalty for :Medium.\n\n\n\n\n\n","category":"type"},{"location":"API/#CellularPotts.VolumePenalty","page":"API","title":"CellularPotts.VolumePenalty","text":"VolumePenalty(λᵥ::Vector{Int})\n\nAn concrete type that penalizes cells that deviate from their desired volume.\n\nRequires a vector λᵥ with n penalties where n is the number of cell types. λᵥ is zero-indexed meaning λᵥ[0] corresponds to the :Medium volume penalty (which is set to zero).\n\nNote: λᵥ is automatically transformed to be a zero-indexed offset array and does not require the volume penalty for :Medium.\n\n\n\n\n\n","category":"type"},{"location":"API/#CellularPotts.addCellProperty-Tuple{CellTable, Symbol, Any, Symbol}","page":"API","title":"CellularPotts.addCellProperty","text":"addCellProperty(df::CellTable, propertyName, defaultValue, cellName)\naddCellProperty(df::CellTable, propertyName, values)\n\nGiven a cellTable, add a new column called propertyName with a given default value.\n\nA cellTable is generated using the newCellState() function.\n\ncellName can be the name of one cell type or a vector of cell types. Cells not included in cellName will be given a value of missing.\n\nIf cellNames are not specified, a vector of values can be supplied for every cell type.  \n\n\n\n\n\n","category":"method"},{"location":"API/#CellularPotts.addNewCell-Union{Tuple{T}, Tuple{CellTable, T}} where T<:NamedTuple","page":"API","title":"CellularPotts.addNewCell","text":"addNewCell(df::CellTable, cell<:NamedTuple)\n\nGiven a cellTable, add a new row corresponding to a new cell in the mode. Property names in the for the cell need to match column names in the cellTable\n\n\n\n\n\n","category":"method"},{"location":"API/#CellularPotts.countcells-Tuple{CellPotts}","page":"API","title":"CellularPotts.countcells","text":"countcells(cpm::CellPotts)\ncountcells(df::CellTable)\n\nCount the number of cells in the model \n\n\n\n\n\n","category":"method"},{"location":"API/#CellularPotts.countcelltypes-Tuple{CellPotts}","page":"API","title":"CellularPotts.countcelltypes","text":"countcelltypes(cpm::CellPotts)\ncountcelltypes(df::CellTable)\n\nCount the number of cell types in the model \n\n\n\n\n\n","category":"method"},{"location":"API/#CellularPotts.newCellState-Union{Tuple{T}, Tuple{Vector{Symbol}, Vector{T}, Vector{T}}} where T<:Integer","page":"API","title":"CellularPotts.newCellState","text":"newCellState(names, volumes, counts)\n\nCreate a new cellTable where each row corresponds to a cell.\n\nBy default, this function generates a table with the following columns:\n\nnames::Vector{Symbol}: List of names given to cells (e.g. :TCell)\ncellIDs::Vector{<:Integer}: A unqiue number given to a cell\ntypeIDs::Vector{<:Integer}: A number corresponding to the cell's name\nvolumes::Vector{<:Integer}: Number of grid squares occupied \ndesiredVolumes::Vector{<:Integer}: Desired number of grid square\nperimeters::Vector{<:Integer}: Cell border penality\ndesiredPerimeters::Vector{<:Integer}: Desired cell border penality\n\nThe first row in the table is reserved for :Medium which is the name given to grid locations not belonging to any cell and is given an index of 0 (The first cell is given an index of 1).\n\nOf note, desiredPerimeters are calculated as the minimal perimeter given the cell's volume. \n\n\n\n\n\n","category":"method"},{"location":"API/#CellularPotts.removeCell-Union{Tuple{T}, Tuple{CellTable, T}} where T<:Integer","page":"API","title":"CellularPotts.removeCell","text":"removeCell(df::CellTable, cellID)\n\nGiven a cellTable, remove the cell with provided cellID.\n\n\n\n\n\n","category":"method"}]
}
