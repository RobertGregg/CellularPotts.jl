using Documenter
using CellularPotts

makedocs(
    sitename = "CellularPotts.jl",
    format = Documenter.HTML(),
    modules = [CellularPotts],
    pages = [
        "Introduction" => "index.md",
        "Examples" => [
            "./ExampleGallery/HelloWorld/HelloWorld.md",
            "./ExampleGallery/LetsGetMoving/LetsGetMoving.md",
            "./ExampleGallery/OnPatrol/OnPatrol.md",
            "./ExampleGallery/LifeAndDeath/LifeAndDeath.md"],
        "API.md"]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/RobertGregg/CellularPotts.jl.git"
)
