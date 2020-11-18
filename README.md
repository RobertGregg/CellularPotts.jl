# Cellular Potts Modeling 

The goal of this project to to create a cellular agent-based model in julia that integrates Potts modeling, Intracellular dynamics with ODEs, and cytokine diffusion with PDEs. The hope is to develop a package that can model immune cell signaling in a multi-scale manner.



Task List

- [x] Get the basic Ising model to work
- [x] Create a function to plot the borders of the cells
  - Maybe look to Agents.jl
  - Could use linesegments + heatmap with makie
- [ ] Introduce cell properties
  - Division
  - Death
  - Active movement 
- [ ] Integrate hybrid modeling schemes
  - [ ] ODE Modeling (intracellular)
  - [ ] PDE Modeling (extracellular)
  - Maybe use [Neural networks](https://github.com/SciML/NeuralPDE.jl) to speed up the PDE computation?
- [ ] Incorporate SDEs for noise in signaling