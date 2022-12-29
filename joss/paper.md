---
title: 'CellularPotts.jl: Simulating Cellular Interactions with Julia'
tags:
  - Julia
  - Cellular Potts
  - Computational Biology
  - Multi-Scale Modeling
authors:
  - name: Robert W. Gregg
    orcid: 0000-0002-2930-621X
    corresponding: true
    affiliation: "1, 2"
  - name: Panayiotis V. Benos
    orcid: 0000-0003-3172-3132
    affiliation: 3
affiliations:
 - name: Department of Computational and Systems Biology, University of Pittsburgh, Pittsburgh, PA, USA
   index: 1
 - name: Department of Biomedical Informatics, University of Pittsburgh, Pittsburgh, PA, USA
   index: 2
 - name: Department of Epidemiology, University of Florida, Gainesville, FL, USA
   index: 3
date: 27 December 2022
bibliography: paper.bib
---

# Summary

CellularPotts.jl is a software package written in Julia to simulate common behaviors seen in biological cells such as division, adhesion, and intercellular signaling. These cellular behaviors influence many complex biological processes (e.g., tumor growth, wound healing, and infection) so if we want to understand these processes, we need software that can accurately model these behaviors. Here we take advantage of Cellular Potts Modeling to simulate cellular behaviors [@graner1992simulation]. These models are advantageous over other approaches because they can simulate changes in cell geometry while still remaining computationally efficient. Users of this package are required to define three key inputs to create a valid model: a 2 or 3 dimensional space, a table describing the cells to be positioned in that space, and a list of model penalties that help dictate cell behavior. Models can then be evolved over time to collect statistics about the observed cell behaviors, simulated repeatedly to understand how a specific property affects the model, and visualized using any of the available plotting libraries in Julia.   

# Statement of need

The goal of this package is to develop Cellular Potts Models (CPMs) in Julia using a graph-based approach. Currently, other excellent software exists (Morpheus [@starruss2014morpheus] Artistoo [@wortel2021artistoo] CompuCell3D [@swat2012multi]) to simulate these types of models, but they have a number of limitations which are addressed by CellularPotts.jl. 

Many of the available packages are written in a low-level language like C++ and are accessed through a GUI or python frontend. This unfortunately separates developers from users, complicates the codebase with multiple languages, and can make user customization difficult. CellularPotts.jl has the advantage of being written purely in Julia: a fast, high-level language tailored to scientific computing. This also comes with the benefit of composing with other Julia packages for easy customization. For example, most CPM software relies on Runge-Kutta or Euler methods to simulate cellular signaling which can easily become unstable and deviate away from the true solution. Julia’s DifferentialEquations.jl library [@rackauckas2017differentialequations] (and more broadly the SciML ecosystem) offers a best-in-class suite of ordinary differential equation (ODE) solvers that can uniquely handle systems where the number of states changes over time and discontinuous jumps occur at any time-point. ODEs can evolve independent of the CPM’s time steps and can even adapt if the error does not meet a specified tolerance.

CellularPotts.jl is distinguished from other CPM software because it defines the space cells occupy as a graph instead of a lattice. Representing models as a graph allows access to decades of graph theory research and solves a number of problems plaguing these types of models. For example, CPMs are notorious for allowing cells to fragment into multiple parts. This is usually addressed by lowering a modeling parameter called "temperature" which only decreases the probability of disconnections from occurring. By encoding cell space as a graph, we can simply test for articulations points [@hopcroft1973algorithm] and avoid locations that would disconnect a cell. This method is guaranteed to work at any model temperature and is independent of how the user defines the topology of the space. 

Another benefit to using a graph-based approach can be seen when modeling cellular division. When a cell divides, the resulting daughter cells are roughly equal in size which is difficult to simulate because cells are irregularly shaped. One method that attempts to solve this issue involves finding a line that optimally divides the cell in half. However, this method is not guaranteed to work if the cell’s boundary is concave. Conversely, by using a graph partitioning algorithm we can handle any cell shape. In general, graph partitioning is an NP hard problem, but because we are only performing a 2-way partition and cells don’t divide with every step, polynomial time heuristic algorithms remain tractable [@karypis1998fast].


# Example

The following code is derived from the "HelloWorld" documentation example ([1]([Hello World · CellularPotts.jl](https://robertgregg.github.io/CellularPotts.jl/dev/ExampleGallery/HelloWorld/HelloWorld/))) which steps through the creation of a model containing one epithelial cell.  

```julia
#Install the package
using Pkg; Pkg.add("CellularPotts")

#Load the package
using CellularPotts

#Create a space (50×50 pixels) for cells to exist in
space = CellSpace(50,50; isPeriodic=true, neighborhood=:moore)

#Describe the cells in the model
initialCellState = CellTable(
    [:Epithelial], #names
    [500],         #sizes
    [1])           #counts

#Add penalties to the model
  #AdhesionPenalty discourages having neighbors with different cell types
  #VolumePenalty discourages cells deviating from their desired size
penalties = [
    AdhesionPenalty([0 20; 20 0]),
    VolumePenalty([5])
    ]

#Create a model object
cpm = CellPotts(space, initialCellState, penalties)

#Record a simulation of the model
recordCPM("Output.gif", cpm)
```

# Acknowledgements

This work was supported by the National Institute of Health (NIH) under the NLM 5T15LM007059-35 grant. 

# References