```@meta
DocTestSetup = quote
    using CellularPotts
end
```

# Introduction

CellularPotts.jl is a package designed to model complex behaviors observed in biological cells. This is accomplished using (you guessed it) Cellular Potts Modeling (CPM). See [What are Cellular Potts Models?](@ref) for more information on how CPMs work. The goal of CellularPotts.jl is to provide a simple interface that is expressive enough to handle diverse behaviors like:

- Adhesion
- Division
- Movement though chemotaxis and contraction of actin filaments
- Cytokine diffusion inside and outside of cells

These behaviors can be combined to develop complex simulations of tumor micro-environments, immune cell activation, stem-cell differentiation, and much more. 

# Installation

```julia
using Pkg
Pkg.add("CellularPotts")
```

# Example Gallery

For a basic introduction, check out the [Hello World](@ref) example or other examples shown below.

!!! tip
    Click on an image to go to example page

```@raw html

<style>
    .Gallery {
        display: grid;
        grid-template-columns: repeat(3,1fr);
        grid-gap: 1%;
        min-width: 0px;
    }  

    h3 {
        text-align: center;
    }
</style>

<div class="Gallery">

    <div>
        <h3>Hello world</h3>
        <a href="https://robertgregg.github.io/CellularPotts.jl/dev/ExampleGallery/HelloWorld/HelloWorld/">
            <img title="" src="https://github.com/RobertGregg/CellularPotts.jl/blob/master/docs/src/ExampleGallery/HelloWorld/HelloWorld.gif?raw=true">
        </a>
    </div>

    <div>
        <h3>Let's get moving</h3>
        <a href="https://robertgregg.github.io/CellularPotts.jl/dev/ExampleGallery/LetsGetMoving/LetsGetMoving/">
            <img title="" src="https://github.com/RobertGregg/CellularPotts.jl/blob/master/docs/src/ExampleGallery/LetsGetMoving/LetsGetMoving.gif?raw=true">
        </a>
    </div>

    <div>
        <h3>On Patrol</h3>
        <a href="https://robertgregg.github.io/CellularPotts.jl/dev/ExampleGallery/OnPatrol/OnPatrol/">
            <img title="" src="https://github.com/RobertGregg/CellularPotts.jl/blob/master/docs/src/ExampleGallery/OnPatrol/OnPatrol.gif?raw=true">
        </a>
    </div>

    <div>
        <h3>Bringing ODEs to life</h3>
        <a href="https://robertgregg.github.io/CellularPotts.jl/dev/ExampleGallery/BringingODEsToLife/BringingODEsToLife/">
            <img title="" src="https://github.com/RobertGregg/CellularPotts.jl/blob/master/docs/src/ExampleGallery/BringingODEsToLife/BringingODEsToLife.gif?raw=true">
        </a>
    </div>

    <div>
        <h3>Going 3D</h3>
        <a href="https://robertgregg.github.io/CellularPotts.jl/dev/ExampleGallery/Going3D/Going3D/">
            <img title="" src="https://github.com/RobertGregg/CellularPotts.jl/blob/master/docs/src/ExampleGallery/Going3D/Going3D.gif?raw=true">
        </a>
    </div>

    <div>
        <h3>Diffusion Outside Cells</h3>
        <a href="https://robertgregg.github.io/CellularPotts.jl/dev/ExampleGallery/DiffusionOutsideCells/DiffusionOutsideCells/">
            <img title="" src="https://github.com/RobertGregg/CellularPotts.jl/blob/master/docs/src/ExampleGallery/DiffusionOutsideCells/DiffusionOutsideCells.gif?raw=true">
        </a>
    </div>

    <div>
        <h3>Diffusion Inside Cells</h3>
        <a href="https://robertgregg.github.io/CellularPotts.jl/dev/ExampleGallery/DiffusionInsideCells/DiffusionInsideCells/">
            <img title="" src="https://github.com/RobertGregg/CellularPotts.jl/blob/master/docs/src/ExampleGallery/DiffusionInsideCells/DiffusionInsideCells.gif?raw=true">
        </a>
    </div>

    <div>
        <h3>Tight Spaces</h3>
        <a href="https://robertgregg.github.io/CellularPotts.jl/dev/ExampleGallery/TightSpaces/TightSpaces/">
            <img title="" src="https://github.com/RobertGregg/CellularPotts.jl/blob/master/docs/src/ExampleGallery/TightSpaces/TightSpaces.gif?raw=true">
        </a>
    </div>

    <div>
        <h3>Over Here</h3>
        <a href="https://robertgregg.github.io/CellularPotts.jl/dev/ExampleGallery/OverHere/OverHere/">
            <img title="" src="https://github.com/RobertGregg/CellularPotts.jl/blob/master/docs/src/ExampleGallery/OverHere/OverHere.gif?raw=true">
        </a>
    </div>

</div>

```