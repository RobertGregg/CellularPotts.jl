# <img title="CellularPotts.jl" src="docs/src/assets/logo.svg" alt="" width="50"> CellularPotts.jl

[docs-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-url]: https://robertgregg.github.io/CellularPotts.jl/dev/

[![][docs-img]][docs-url] [![codecov](https://codecov.io/gh/RobertGregg/CellularPotts.jl/graph/badge.svg?token=D3GKFH900T)](https://codecov.io/gh/RobertGregg/CellularPotts.jl) [![CI](https://github.com/RobertGregg/CellularPotts.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/RobertGregg/CellularPotts.jl/actions/workflows/CI.yml)

**CellularPotts.jl** is a Julia package designed to simulate behaviors observed in biological cells like division and adhesion. Users of this package can create 2D and 3D environments with any number of cell types, sizes, and behaviors. Simulations can be recorded and visualized as animations with the help of the Plots.jl package. The goal of this package is to create a flexible coding environment to explore how cell behaviors can coalesce into complex dynamics while still maintaining high performance. Compared to other excellent software for Cellular Potts modeling (e.g., [Morpheus](https://morpheus.gitlab.io/), [Artistoo](https://artistoo.net/), [CompuCell3D](https://compucell3d.org/)), CellularPotts.jl is unique in its approach for a few reasons:

- CellularPotts.jl is written completely in Julia, avoiding the "[two language problem](https://www.nature.com/articles/d41586-019-02310-3)"
  
  - This unites developers and users to one language, simplifies the code base, and makes customization easier.

- The space cells occupy is modeled as a network/graph
  
  - Representing the model as a graph allows access to decades of graph theory research, for example:
    
    - Calculating articulation points to avoid cell fragmentation
    - Using graph partitioning algorithms to simulate cellular division
    - Avoiding cumbersome boundary conditions by adding edges that wrap around the defined space
    - Using graphical Laplacians to simulate diffusion

- CellularPotts.jl can be composed with other Julia packages.
  
  - For example, we can use state-of-the-art differential equation solving techniques from [DifferentialEquations.jl](https://diffeq.sciml.ai/stable/) as opposed to simple Euler methods. (Check out the [BringingODEStoLife]([Bringing ODEs To Life · CellularPotts.jl](https://robertgregg.github.io/CellularPotts.jl/dev/ExampleGallery/BringingODEsToLife/BringingODEsToLife/)) example.)

## Quick Start

To create a basic Cellular Potts Model, you need to provide 3 pieces of information:

1. What space will the cells occupy?

2. What cells do you want to include in the model?

3. What penalties do you want to encourage certain behaviors?

```julia
#Install the package (if needed)
using Pkg; Pkg.add("CellularPotts")

#Load in the package
using CellularPotts

#Create a space (50×50) for cells to exist in
space = CellSpace(50,50; periodic=true, diagonal=true)

#Describe the cells in the model
initialCellState = CellState(
    :Epithelial,         #cell names
    500,                 #size of cells
    1,                   #number of cells
    positions = (25,25))

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

## Want to Contribute?

Careful attention has been taken to ensure this package is as performant as I can possibly make it, however, if you spot something egregious in the package, feel free to raise an issue or pull request.

Also of note, **this package is still in development and is not currently recommended for general use**. However, still feel free to try it and give suggestions if you're curious. 

## How to Cite

If you happen to use CellularPotts.jl in your research, feel free to cite our paper:

https://doi.org/10.1093/bioinformatics/btad773

```
@article{gregg2024cellularpotts,
  title={CellularPotts. jl: simulating multiscale cellular models in Julia},
  author={Gregg, Robert W and Benos, Panayiotis V},
  journal={Bioinformatics},
  volume={40},
  number={1},
  pages={btad773},
  year={2024},
  publisher={Oxford University Press}
}
```

## Future Improvements

- [ ] `CellDivision!()` currently cannot update custom cell state properties

- [ ] Keyword options for cell state (to add cell properties)

- [ ] Use automatic differentiation to calculate [cellular forces](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1007459) from the Hamiltonian

- [x] Create more unit tests for reproducibility (low code coverage)

- [ ] Use SVectors to store graph edges? 🤔
  
  - Only useful for spaces where all nodes are identical (e.g., periodic boundaries)

- [x] Add CI and CodeCov badge

- [x] Reduce heavy package dependencies using package extensions
  
  - [ ] Not using package extensions yet but some heavy dependencies were removed

- [ ] Separate tutorials from examples

- [ ] Tutorial on how to create your own penalty
