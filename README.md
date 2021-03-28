# Cellular Potts Modeling 

The goal of this project to to create a cellular agent-based model in julia that integrates Potts modeling, Intracellular dynamics with ODEs, and cytokine diffusion with PDEs. The hope is to develop a package that can model immune cell signaling in a multi-scale manner.



## Task List

- [x] Get the basic Ising model to work
- [x] Create a function to plot the borders of the cells
  - Could use linesegments + heatmap with makie
- [ ] Introduce cell properties
  - [x] Division
  - [ ] Death
  - [ ] Active movement
  - [ ] Movement up gradients
- [ ] Integrate hybrid modeling schemes
  - [ ] ODE Modeling (intracellular)
  - [ ] PDE Modeling (extracellular)
  - Maybe use [Neural networks](https://github.com/SciML/NeuralPDE.jl) to speed up the PDE computation?
- [ ] Incorporate SDEs for noise in signaling
- [ ] Hook into [Agents.jl](https://github.com/JuliaDynamics/Agents.jl) or [DynamicGrids.jl](https://github.com/cesaraustralia/DynamicGrids.jl) ?

## For cell gradient:

- visualization of the gradient, can you overlay heatmaps (top one being transparent?)
  - you can have overlapping heatmap! see [deletethis.jl](./src/deletethis.jl)
- How do you store the gradient field? as another field in CellPotts?
- What if you need more than one gradient field?
- How will the gradient act "inside" a cell? Will it be constant there?
- How will cells update the gradient? Just at their border?

## For cell division

  - Should the cell age be tracked just like the volume?
  - ~~Use the center/line split method or the graph partition method?~~ 
      - decided to use center/line split method, optimizing the slope to most evenly divide the cell
  - ~~Need to update a lot of variables after the division (grid, volume, lambda)~~ 
      - just need to make the H re-calculation its own function

## TODO Fixes

- [ ] If the grid wraps around, the Edge2Grid function will have repeated edges
- [ ] separate a struct for a cell and for a grid? Would having a grid of cells be slower/more memory intensive?
- [ ] cellDivideButton seems to divide the cells twice before the GUI starts