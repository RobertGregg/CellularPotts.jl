```@meta
EditURL = "<unknown>/docs/src/ExampleGallery/Going3D/Going3D.jl"
```

# Going 3D

````julia
using CellularPotts


space = CellSpace(30,30,30; wrapAround=false)


initialCellState = CellTable(
    [:Epithelial],
    [1000],
    [1])


positions = [size(space).÷2]

initialCellState = addcellproperty(initialCellState, :positions, positions)

penalties = [
    AdhesionPenalty([0 30;
                     30 0]),
    VolumePenalty([5])
    ]


cpm = CellPotts(space, initialCellState, penalties)

positionCells!(cpm)

recordCPM("Going3D.gif", cpm)
````

```@raw html
<p style="text-align:center;">
    <img title="" src="https://github.com/RobertGregg/CellularPotts.jl/blob/master/docs/src/ExampleGallery/Going3D/Going3D.gif?raw=true" width="445">
</p>
```


---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

