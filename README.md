# NEED TO FIX CELL POTTS MODEL 
- choose a random lattice site i
- choose a random neighboring lattice site j to copy its ID into i. (THIS NEEDS TO BE FIXED)

# Cellular Potts Modeling 

The goal of this project to to create a cellular agent-based model in julia that integrates Potts modeling, Intracellular dynamics with ODEs, and cytokine diffusion with PDEs. The hope is to develop a package that can model immune cell signaling in a multi-scale manner.



Task List

- [x] Get the basic Ising model to work
- [x] Create a function to plot the borders of the cells
  - Maybe look to Agents.jl
  - Could use linesegments + heatmap with makie
- [ ] Introduce cell properties
  - [ ] Division
  - [ ] Death
  - [ ] Active movement
  - [ ] Movement up gradients
- [ ] Integrate hybrid modeling schemes
  - [ ] ODE Modeling (intracellular)
  - [ ] PDE Modeling (extracellular)
  - Maybe use [Neural networks](https://github.com/SciML/NeuralPDE.jl) to speed up the PDE computation?
- [ ] Incorporate SDEs for noise in signaling


For cell gradient:
- visualization of the gradient, can you overlay heatmaps (top one being tranparent?)
- How do you store the gradient field? as another field in CellPotts?
- What if you need more than one gradient field?
- How will the gradient act "inside" a cell? Will it be contant there?
- How will cells update the gradient? Just at their border?

For cell division
  - Should the cell age be tracked just like the volume?
  - Use the center/line split method or the graph partition method?
  - Need to update a lot of variables (grid, volume, lambda) 