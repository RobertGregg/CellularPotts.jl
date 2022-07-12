using CellularPotts
using GLMakie


M = ModelParameters(
    graphDimension = (100,100),
    isPeriodic = true,
    cellTypes = ["SimpleCell"],
    cellCounts = [1],
    cellVolumes = [500],
    penalties = [
        AdhesionPenalty([0 20; 20 50]),
        VolumePenalty([500],[1])
    ],
    temperature = 10.0
)

CPM = CellPotts(M)

CellGUI(CPM)