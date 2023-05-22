# # Over Here!
using Revise
using CellularPotts

space = CellSpace(50,50)


initialCellState = CellTable(
    [:Epithelial],
    [200],
    [1])


positions = [size(space).÷2]
    
initialCellState = addcellproperty(initialCellState, :positions, positions)

#Species concentration gradient
species = [100exp(-(x^2+y^2)/1000) for x in first(axes(space)), y in last(axes(space))]

penalties = [
    AdhesionPenalty([0 30;
    30 0]),
    VolumePenalty([5]),
    ChemoTaxisPenalty([1], species)
]


cpm = CellPotts(space, initialCellState, penalties)



recordCPM("OverHere.gif", cpm, timestamps=0:3000, framerate=60)

