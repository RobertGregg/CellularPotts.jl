```@meta
EditURL = "Going3D.jl"
```

# Going 3D

````julia
using CellularPotts


space = CellSpace(30,30,30; periodic=false)


initialCellState = CellState(names=:Epithelial, volumes=1000, counts=2, positions = [(12,12,12),(18,18,18)])

penalties = [
    AdhesionPenalty([0 20;
                     20 20]),
    VolumePenalty([1])
    ]


cpm = CellPotts(space, initialCellState, penalties)

record(cpm; file="Going3D.gif", colorby=:id, cellcolors = [:dodgerblue, :red], size=(800,800))
````

```@raw html
<p style="text-align:center;">
    <img title="" src="https://github.com/RobertGregg/CellularPotts.jl/blob/master/docs/src/ExampleGallery/Going3D/Going3D.gif?raw=true" width="445">
</p>
```


---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

