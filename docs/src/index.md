```@meta
DocTestSetup = quote
    using CellularPotts
end
```

# Introduction

CPMs simulate a collection of cells interacting with one another. These interactions can range anywhere from simple contact between cells to complex cytokine communication.

# Installation

```julia
using Pkg
Pkg.add("CellularPotts")
```

# Example Gallery

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

</div>

```