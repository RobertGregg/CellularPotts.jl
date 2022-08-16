using Documenter
using CellularPotts


makedocs(
    sitename = "CellularPotts.jl",
    format = Documenter.HTML(),
    modules = [CellularPotts],
    pages = [
        "Introduction" => "index.md",
        "Examples" => [
            "Hello World" => "ExampleGallery/HelloWorld/HelloWorld.md",
            "Let's Get Moving" => "ExampleGallery/LetsGetMoving/LetsGetMoving.md",
            "On Patrol" => "ExampleGallery/OnPatrol/OnPatrol.md",
            "Bringing ODEs To Life" => "ExampleGallery/BringingODEsToLife/BringingODEsToLife.md",
            "Going 3D" => "ExampleGallery/Going3D/Going3D.md"],
        "API.md"]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/RobertGregg/CellularPotts.jl.git"
)
