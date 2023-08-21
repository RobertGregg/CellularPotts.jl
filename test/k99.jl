using CellularPotts

#Diameter of airway: 1000 μm
airwayDiameter = 1000
airwayRadius = airwayDiameter/2
airwayArea = π*airwayRadius^2 |> x->round(Int,x)

#Diamter of respiratory epithelial cell: 20 μm
cellDiameter = 20
cellArea = π*(cellDiameter/2)^2 |> x->round(Int,x)
totalCells = 700

boxSize = 1200

#Make a space with hole in middle
spaceImage = zeros(Int, boxSize,boxSize)
for i in axes(spaceImage,2)
    for j in axes(spaceImage,1)
        if airwayRadius^2 ≤ (i-boxSize ÷ 2)^2 + (j-boxSize ÷ 2)^2 ≤ (airwayRadius+59)^2
            spaceImage[i,j] = 1
        end
    end
end

space = CellSpace(spaceImage, periodic=false)

initialCellState = CellState(
    [:Ciliated, :Club, :Basal],
    [cellArea,cellArea,cellArea÷2],
    [500,100,100]
    );

airwayPosition = size(space) .÷ 2
positions  = Vector{Tuple{Int64, Int64}}()
cellsPositioned = 0

while cellsPositioned ≤ totalCells
    r = sqrt(rand()) * boxSize/2
    θ = 2π*rand()

    if airwayRadius ≤ r ≤ airwayRadius+5
        cellPosition = (r*cos(θ), r*sin(θ)) .+ airwayPosition .|> x->round(Int,x)
        push!(positions, cellPosition)
        cellsPositioned += 1
    end
end



initialCellState = addcellproperty(initialCellState, :positions, positions)

penalties = [
    AdhesionPenalty([20 30 30 30;
                     30 20 30 30;
                     30 30 20 30;
                     30 30 30 20]),
    VolumePenalty([30, 30, 30, 30]),
    ]


cpm = CellPotts(space, initialCellState, penalties)


# for i in 1:100
#     println(i)
#     ModelStep!(cpm)
# end

p = plotcpm(cpm,property=:nodeTypes)

savefig(p,"airway.pdf", height=3,width=3)


(rows,columns) = size(cpm.space)

plotObject = heatmap(
            getproperty(cpm.space, :nodeTypes)',
            c = cgrad(:tol_light, rev=true),
            grid=false,
            axis=nothing,
            legend=:none,
            framestyle=:box,
            aspect_ratio=:equal,
            size = (600,600),
            xlims=(0.5, rows+0.5),
            ylims=(0.5, columns+0.5),
            clim=(0,3),
            dpi=600
            )

cellborders!(plotObject, cpm.space)

savefig(plotObject,"airway.png")