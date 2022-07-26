using Documenter
using CellularPotts

makedocs(
    sitename = "CellularPotts",
    format = Documenter.HTML(),
    modules = [CellularPotts]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/RobertGregg/CellularPotts.jl.git"
)
