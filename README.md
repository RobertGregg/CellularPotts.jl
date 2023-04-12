# <img title="CellularPotts.jl" src="docs/src/assets/logo.svg" alt="" width="50"> CellularPotts.jl

[docs-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-url]: https://robertgregg.github.io/CellularPotts.jl/dev/

[![][docs-img]][docs-url]

The goal of this package is to develop a Cellular Potts modeling approach in Julia using a network-based approach. Currently, other software exists to simulate these types of models, but they have a number of limitations:

- They are written in a low-level language (e.g. C++) with a GUI or python frontend
  - This separates developers from users, complicates the code base, and makes customization difficult
- They rely on a grid approach instead of a network based approach
  - Representing the model as a graph allows access to decades of graph theory research, for example:
    - calculating articulation points to avoid cells disconnecting
    - using graph partitioning algorithms to simulate cellular division
    - avoiding cumbersome boundary conditions by simply adding edges that loop around
    - using graphical laplacians to simulate diffusion 
- They cannot take advantage of how composable Julia packages are with one another. For example, we can use state-of-the-art differential equation solving techniques from [DifferentialEquations.jl](https://diffeq.sciml.ai/stable/). 
  - Most CPM software relies on Runge-Kutta or even simple Euler Methods

Researchers and developers have been able to accomplish a lot with their respective softwares and I would urge anyone to check them out. My favorites are [Morpheus](https://morpheus.gitlab.io/), [Artistoo](https://artistoo.net/), and [CompuCell3D](https://compucell3d.org/). This package takes a lot of inspiration from their design and pedagogical examples.

Careful attention has been taken to ensure this package is as performant as I can possibly make it, however, if you spot something egregious in the package, feel free to raise an issue or pull request.

Also of note, **this package is still in major development and is not currently recommended for general use**. I'm still working out how to best organize datastructures and functionally. However, still feel free to try it if you're curious. 

## Simple Example

```julia
#Install the package
using Pkg; Pkg.add("CellularPotts")

#Load in the package
using CellularPotts

#Create a space (50×50) for cells to exist in
space = CellSpace(50,50; isPeriodic=true, neighborhood=:moore)

#Describe the cells in the model
initialCellState = CellTable(
    [:Epithelial], #names
    [500],         #sizes
    [1])           #counts

#Add penalties to the model
penalties = [
    AdhesionPenalty([0 20;
                     20 0]),
    VolumePenalty([5])
    ]

#Create a model object
cpm = CellPotts(space, initialCellState, penalties)

#Record a simulation of the model
recordCPM("ReadMeExample.gif", cpm)
```

<img title="ReadMeEaxmple" src="https://github.com/RobertGregg/CellularPotts.jl/blob/master/docs/src/ExampleGallery/HelloWorld/HelloWorld.gif?raw=true">

## Major Improvements

- [ ] Introduce more cell properties
  
  - [x] Division
  - [x] Death
  - [x] Active movement
  - [ ] Movement up gradients

- [x] Integrate hybrid modeling schemes
  
  - [x] ODE Modeling (intracellular)
  
  - [x] PDE Modeling (extracellular)
  
  - Maybe use [Neural networks](https://github.com/SciML/NeuralPDE.jl) to speed up the PDE computation?
  
  - Stochastic jumps?

- [x] Create an Examples folder

- [x] How to save output?
  
  - Save the data into a dictionary of dataframes
  - Needs to be made more efficient

- [x] Implement different ways to initialize cell locations
  
  - [x] Image input
  - [x] specify locations with property

- [x] Allow cells to have different properties (used `NamedTuple`)

- [ ] Use automatic differentiation to calculate [cellular forces](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1007459) from the Hamiltonian

- [x] Add a correction factor to adhesion to deal with boundaries.

## Minor Improvements

- [x] Allow user defined parameters to cells (used `NamedTuple`)
- [x] Allow cells of the same type to be different sizes (?)
  - Just specify different desired volumes
- [x] Could get a big speed improvement if you don't loop through all cells to update articulation points
  - Need to be clever about updating articulation points locally (is this possible?)
  - rewrote [Tarjan's algoirthm](https://en.wikipedia.org/wiki/Biconnected_component) to find articulation points which is O(V+E)
- [x] Adding cell borders is slow for large spaces
  - fixed by using NA
- [ ] Use abstract typing (e.g. `AbstractVector` vs `Vector`) without creating type instability
- [x] Do we even need to track to total energy? (nope!)
- [ ] Use SVectors to store graph edges? 🤔
  - Only useful for spaces where all nodes are identical (e.g., periodic boundaries)
- [ ] Add more tests and CI badge