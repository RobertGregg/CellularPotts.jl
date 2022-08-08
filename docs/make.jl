using Documenter
using CellularPotts

makedocs(
    sitename = "CellularPotts.jl",
    format = Documenter.HTML(),
    modules = [CellularPotts],
    pages = [
        "Introduction" => "index.md",
        "Examples" => ["./GeneratedGallery/HelloWorld.md"],
        "API.md"]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/RobertGregg/CellularPotts.jl.git"
)
