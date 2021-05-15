# Cellular Potts Modeling 

The goal of this project to to create a cellular agent-based model in julia that integrates Potts modeling, Intracellular dynamics with ODEs, and cytokine diffusion with PDEs. The hope is to develop a package that can model immune cell signaling in a multi-scale manner.



## Questions

- Best way to implement gradient field?
  - How will the gradient act "inside" a cell? Will it be constant there?
  - How will cells update the gradient? Just at their border?
- Should the cell age be tracked for division?
- What is the difference b/w Rect3D and FRect3D?

## Comments

- An ⊗ I(m) + I(n) ⊗ Am ≠ Am ⊗ I(n) + I(m) ⊗ An when n≠m
- [Metagraphs.jl](https://github.com/JuliaGraphs/MetaGraphs.jl) saves attributes as `Dict{Symbol, Any}` which leads to a lot of type instability
  - The upside is you can dynamically add any number of attributes
- Cells on opposite borders without periodic boundaries will disconnect medium
- If you see an offset array, its just to give medium a zero index (so cell 1 can have index 1)

## Major Improvements

- [ ] Introduce cell properties
  - [x] Division
  - [ ] Death
  - [ ] Active movement
  - [ ] Movement up gradients
- [ ] Integrate hybrid modeling schemes
  - [ ] ODE Modeling (intracellular)
  - [ ] PDE Modeling (extracellular)
  - Maybe use [Neural networks](https://github.com/SciML/NeuralPDE.jl) to speed up the PDE computation?
- [ ] Hook into [Agents.jl](https://github.com/JuliaDynamics/Agents.jl) or [DynamicGrids.jl](https://github.com/cesaraustralia/DynamicGrids.jl) ?
- [ ] Create an Examples folder
- [ ] How to save output?

## Minor Improvements

- [ ] Allow additional parameters to cells (maybe input as a named tuple?)
- [ ] Allow cells of the same type to be different sizes
- [ ] For 3D gui, don't recreate voxels every iteration, just set color to clear
- [ ] Could get a big speed improvement if you don't loop through all cells to update articulation points
  - need to be clever about updating articulation points locally
- [ ] VP having desired volumes makes cell division difficult (type unstable)
  - replace with CPM.M.cellVolumes ?