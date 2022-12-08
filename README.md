
[docs-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-url]: https://robertgregg.github.io/CellularPotts.jl/dev/

[![][docs-img]][docs-url]

# <img title="CellularPotts.jl" src="docs/src/assets/logo.svg" alt="" width="50"> CellularPotts.jl

The goal of this package is to develop a cellular, agent-based modeling approach in Julia using a network-based Potts modeling framework. Currently, other software exists to simulate these types of models, but they have a number of limitations:

- They are written in a low-level language (e.g. C++) with a GUI or python frontend
  - This separates developers from users, complicates the code base, and makes customization difficult
- They rely on a grid approach instead of a network based approach
  - representing the model as a graph allows access to decades of graph theory research, for example:
    - calculating articulation points to avoid cells disconnecting
    - using graph partitioning algorithms to simulate cellular division
    - avoiding cumbersome boundary conditions by simply adding edges that loop around
- They cannot take advantage of how composable Julia packages are with one another. For example, we can use state-of-the-art differential equation solving techniques from [DifferentialEquations.jl](https://diffeq.sciml.ai/stable/). 
  - Most CPM software relies on Runge-Kutta or even simple Euler Methods

Researchers and developers have been able to accomplish a lot with their respective softwares and I would urge anyone to check them out. My favorites are [Morpheus](https://morpheus.gitlab.io/), [Artistoo](https://artistoo.net/), and [CompuCell3D](https://compucell3d.org/). This package takes a lot of inspiration from their design and pedagogical examples.

Careful attention has been taken to ensure this package is as performant as I can possibly make it (particularly with type stability and allocations). For example, stepping the model forward in time produces no allocations. However, if you spot something egregious in the codebase, feel free to raise an issue or pull request.

Also of note, **this package is still in major development and is not currently recommended for general use**. I'm still working out how to best organize datastructures and functionally. However, still feel free to try it if you're curious. 

Below is simply tracking the progress of the package and any notes to myself.

## Questions

- Best way to implement gradient field?
  - How will the gradient act "inside" a cell? Will it be constant there?
  - How will cells update the gradient? Just at their border?
- Should cell age be tracked for division?
- What is the difference b/w Rect3D and FRect3D?
- How to hook into Agents.jl?

## Comments

- [Metagraphs.jl](https://github.com/JuliaGraphs/MetaGraphs.jl) saves attributes as `Dict{Symbol, Any}` which leads to a lot of type instability
  - The upside is you can dynamically add any number of attributes
- A wall of cells expanding across a grid without periodic boundaries will disconnect medium (location with no cell present)
- Offset arrays are used to give Medium a zero index (so cell 1 can have index 1). I'm not use if this will lead to issues later on or leads to a loss in proformance. 

## Major Improvements

- [ ] Introduce more cell properties
  
  - [x] Division
  - [x] Death
  - [x] Active movement
  - [ ] Movement up gradients

- [ ] Integrate hybrid modeling schemes
  
  - [x] ODE Modeling (intracellular)
  
  - [ ] PDE Modeling (extracellular)
  
  - Maybe use [Neural networks](https://github.com/SciML/NeuralPDE.jl) to speed up the PDE computation?
  
  - Stochastic jumps?

- [x] Create an Examples folder

- [x] How to save output?
  - Save the data into a dictionary of dataframes
  - Needs to be made more efficient

- [x] Implement different ways to initialize cell locations
  
  - [ ] Still need image input

- [x] Allow cells to have different properties (used `NamedTuple`)

- [ ] Use automatic differentiation to calculate [cellular forces](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1007459) from the Hamiltonian

## Minor Improvements

- [x] Allow user defined parameters to cells (used `NamedTuple`)
- [x] Allow cells of the same type to be different sizes (?)
  - Just specify different desired volumes
- [ ] For 3D gui, don't recreate voxels every iteration, just set color to clear
- [x] Could get a big speed improvement if you don't loop through all cells to update articulation points
  - Need to be clever about updating articulation points locally (is this possible?)
  - rewrote [Tarjan's algoirthm](https://en.wikipedia.org/wiki/Biconnected_component) to find articulation points which is O(V+E)
- [ ] Don't loop through all of the cell borders when updating the GUI.
- [ ] Use abstract typing (e.g. `AbstractVector` vs `Vector`) without creating type instability
- [x] Do we even need to track to total energy? (nope!)
