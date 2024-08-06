using CellularPotts
using GLMakie
using Colors


mutable struct Visual
    figure::Figure
    axis::Axis
    hm::Heatmap{Tuple{Vector{Float32}, Vector{Float32}, Matrix{Float32}}}
    ls::Lines{Tuple{Vector{Point{2, Float32}}}}
    sizeSpace::Tuple{Int64,Int64}
    cellColors::Vector{RGBA{Float64}}
    migrationColors::Vector{RGBA{Float64}}
    borders::Vector{Point{2, Float32}}

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
        migrationColors = range(RGBA(255,255,255,0.0), RGBA(60,4,104,0.6))
    
        sizeSpace = size(cpm.space)
        hm = heatmap!(axis, cpm.space.nodeIDs; colormap=Makie.Categorical(cellColors))

        borders = Point2f.(fill(NaN, 3*prod(sizeSpace)),fill(NaN, 3*prod(sizeSpace)))
        updateBorders!(borders, cpm.space)

        ls = lines!(axis, borders, color=:black)

    
        visual =  new(
            figure,
            axis,
            hm,
            ls,
            sizeSpace,
            cellColors,
            migrationColors,
            borders
        )


        return visual
    end
end

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



space = CellSpace(50,50; periodic=false, diagonal=true)
initialCellState = CellState([:A,:B,:C,:D], [425,425,425,425], [1,1,1,1]; positions = [(10,10),(40,40),(10,40),(40,10)])
penalties = [
    AdhesionPenalty([0 20 20 20 20;
                     20 0 20 20 20;
                     20 20 0 20 20;
                     20 20 20 0 20;
                     20 20 20 20 0]),
    VolumePenalty([5,5,5,5])
    ]

cpm = CellPotts(space, initialCellState, penalties)

visual = Visual(cpm)

record(visual.figure, "test.gif", 1:300) do i
    ModelStep!(cpm)


    updateBorders!(visual, cpm.space)
    
    visual.hm[1] = cpm.space.nodeIDs
    visual.ls[1] = visual.borders
end