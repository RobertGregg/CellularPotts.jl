using CellularPotts
using GLMakie
using Colors

mutable struct Visual
    figure::Figure
    axis::Axis
    sizeSpace::Tuple{Int64,Int64}
    cellColors::Vector{RGBA{Float64}}
    migrationColors::Vector{RGBA{Float64}}
    borders::Vector{Point{2, Float32}}
    migrationIndex::Int64
end


function Visual(cpm::CellPotts)
       
    figure  = Figure(size = (600, 600))

    axis = Axis(figure[1, 1], aspect = DataAspect())
    hidedecorations!(axis)

    cellColors = [
        RGBA(221/255, 221/255, 221/255, 1.0) #off-white
        RGBA(119/255, 170/255, 221/255, 1.0) #light blue
        RGBA(238/255, 221/255, 136/255, 1.0) #pale yellow
        RGBA(68/255, 187/255, 153/255, 1.0) #mint green
        RGBA(238/255, 136/255, 102/255, 1.0) #peach
        RGBA(255/255, 170/255, 187/255, 1.0) #pink
        ]

    cellColors = cellColors[1:(countcelltypes(cpm)+1)]
    migrationColors = range(RGBA(1,1,1,0.0), RGBA(60/255, 4/255, 104/255, 0.6))

    #Cell IDs
    sizeSpace = size(cpm.space)
    heatmap!(axis, cpm.space.nodeIDs; colormap=Makie.Categorical(cellColors))

    #Cell Borders
    borders = Point2f.(fill(NaN, 3*prod(sizeSpace)),fill(NaN, 3*prod(sizeSpace)))
    updateBorders!(borders, cpm.space)

    lines!(axis, borders, color=:black)

    #Cell Migration
    migrationIndex = findfirst(x->x isa MigrationPenalty, cpm.penalties)

    if !isnothing(migrationIndex)
        heatmap!(axis, cpm.penalties[migrationIndex].nodeMemory; colormap=migrationColors)
    else
        migrationIndex = 0
    end

    return Visual(
        figure,
        axis,
        sizeSpace,
        cellColors,
        migrationColors,
        borders,
        migrationIndex
    )

end

Base.getindex(visual::Visual,i) = visual.axis.scene.plots[i]
Base.length(visual::Visual) = length(visual.axis.scene.plots)

updateBorders!(visual::Visual,space::CellSpace) = updateBorders!(visual.borders, space)

function updateBorders!(borders,space::CellSpace)

    (row,col) = size(space)

    
    borderEnd = 0
    for r in 1:row
        for c in 1:col
                        
            #vertical
            if space.nodeIDs[r,c] ≠ space.nodeIDs[r,mod1(c+1,col)]

                borders[borderEnd+1] = Point2f(r-0.5, c+0.5)
                borders[borderEnd+2] = Point2f(r+0.5, c+0.5)

                borderEnd += 3
                
            end

            #horizontal
            if space.nodeIDs[r,c] ≠ space.nodeIDs[mod1(r+1,row),c]

                borders[borderEnd+1] = Point2f(r+0.5, c-0.5)
                borders[borderEnd+2] = Point2f(r+0.5, c+0.5)

                borderEnd += 3
                
            end

        end
    end

    for i in (borderEnd+1):length(borders)
        borders[i] = Point2f(NaN, NaN)
    end

    return nothing
end

space = CellSpace(100,100; diagonal=true)
initialCellState = CellState(:Epithelial, 500, 1; positions = size(space) .÷ 2);


penalties = [
    AdhesionPenalty([0 30;
                    30 30]),
    VolumePenalty([5]),
    PerimeterPenalty([0,10]),
    MigrationPenalty(50, [50], size(space))
    ]


cpm = CellPotts(space, initialCellState, penalties)

visual = Visual(cpm)

record(visual.figure, "test.gif", 1:200) do i

    ModelStep!(cpm)

    updateBorders!(visual, cpm.space)
    
    visual[1][1] = cpm.space.nodeIDs
    visual[2][1] = visual.borders

    if visual.migrationIndex > 0
        visual[3][1] = cpm.penalties[visual.migrationIndex].nodeMemory
    end
end