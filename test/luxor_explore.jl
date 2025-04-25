using Luxor
using CellularPotts

Base.mod1(I::CartesianIndex{N}, i::NTuple{N, Int64}) where N = CartesianIndex(mod1.(Tuple(I),i))


function addAlpha(colormap; minAlpha=0.0, maxAlpha=0.8)
    
    cmap = colorschemes[colormap]
    alphas = range(minAlpha, maxAlpha, length=length(cmap))

    return ColorScheme([RGBA(c.r, c.g, c.b, alpha) for (c,alpha) in zip(cmap,alphas)])
end



function drawHeatmap(Xt; cellsize=(10,10), colormap=:thermal, alpha=0.5)
    
    X = Xt'
    table = Table(size(X)..., cellsize...)

    if !isone(alpha)
        matrixColors = addAlpha(colormap)
    end

    xmin, xmax = extrema(X)
    xrange = xmax - xmin

    for i in LinearIndices(X)
        colorValue = (X[i] - xmin) / xrange

        setcolor(get(matrixColors, colorValue))
        box(table,i, action=:fill)
    end
end


function drawCells(cpm::CellPotts{T,C,2,S,U};
    colorby=:type,
    colormap=nothing,
    migrationcolors=range(colorant"#ffffff00", colorant"#3c0468cc"),
    cellsize=(10,10)) where {T,C,S,U}
    
    #Pick a default color map depending on inputs
    if isnothing(colormap)
        if colorby == :type && countcelltypes(cpm) ≤ 6
            cellcolors = colorschemes[:tol_light]
        else
            cellcolors = get(colorschemes[:darkrainbow], range(0, 1, length=countcells(cpm)))
        end
    else
        cellcolors = get(colorschemes[colormap], range(0, 1, length=countcells(cpm)))
    end

    #Choose a matrix for coloring and its size
    colorMatrix = colorby == :type ? cpm.space.nodeTypes : cpm.space.nodeIDs
    nodeIDs = cpm.space.nodeIDs
    gridSize = size(cpm.space)

    #Create a Luxor table to loop through
    table = Table(gridSize..., cellsize...)
    box(BoundingBox(table), action = :stroke)
    
    up = CartesianIndex(1,0)
    down = CartesianIndex(-1,0)
    left = CartesianIndex(0,-1)
    right = CartesianIndex(0,1)

    for I in CartesianIndices(colorMatrix)
        r, c = Tuple(I)

        #Don't color empty space
        if iszero(colorMatrix[I])
            continue
        end

        b1 = box(table, r, c)


        if colorby ≠ :none
            sethue(cellcolors[colorMatrix[I]])
            poly(b1, action=:fill)
        end

        #Border lines
        sethue("black")
        if nodeIDs[I] ≠ nodeIDs[mod1(I+up, gridSize)]
            line(b1[1],b1[4], action=:stroke)
        end

        if nodeIDs[I] ≠ nodeIDs[mod1(I+left, gridSize)]
            line(b1[1],b1[2], action=:stroke)
        end

        if nodeIDs[I] ≠ nodeIDs[mod1(I+down, gridSize)]
            line(b1[2],b1[3], action=:stroke)
        end

        if nodeIDs[I] ≠ nodeIDs[mod1(I+right, gridSize)]
            line(b1[3],b1[4], action=:stroke)
        end
    end
end




cpm = CellPotts(
        CellSpace(100,100, periodic=true, diagonal=false),
        CellState(names=[:A,:B,:C], volumes=[100,100,100], counts=[10,10,10]),
        [AdhesionPenalty(fill(10,4,4)), VolumePenalty([2,2,2])]
    )

ModelStep!(cpm)

@drawsvg begin
    background("grey90")
    setline(1.0)
    drawCells(cpm; colorby=:type, cellsize = (5,5))
end


X = sort(rand(100,100), dims=2)
@drawsvg begin
    background("grey90")
    sethue("steelblue")
    circle(O,100, action=:fill)
    drawHeatmap(X; colormap=:Purples, cellsize=(5,5))
end





space = CellSpace(100,100; diagonal=true)
initialCellState = CellState(names=:Epithelial, volumes=500, counts=1, positions=size(space) .÷ 2);
penalties = [
    AdhesionPenalty([0 30;
                    30 30]),
    VolumePenalty([5]),
    PerimeterPenalty([0,10]),
    MigrationPenalty(50, [50], size(space))
    ]

cpm = CellPotts(space, initialCellState, penalties)

for i=1:100
    ModelStep!(cpm)
end


migrationIndex = findfirst(x->x isa MigrationPenalty, cpm.penalties)
nodeMemory = cpm.penalties[migrationIndex].nodeMemory
@drawsvg begin
    background("grey90")
    setline(1.0)
    drawCells(cpm; colorby=:type, cellsize = (5,5))
    drawHeatmap(nodeMemory; colormap=:Purples, cellsize=(5,5))
end





mymovie = Movie(400, 400, "mymovie")

function frame(scene::Scene, framenumber::Int64)
    background("white")
    norm_framenumber = rescale(framenumber,
        scene.framerange.start,
        scene.framerange.stop,
        0, 1)
    rotate(norm_framenumber * 2π)
    juliacircles(100)
end

animate(mymovie,
        [
            Scene(mymovie, frame, 1:60)
        ],
    creategif=true,
    pathname="juliaspinner.gif")




mymovie = Movie(400, 400, "mymovie", 1:300)

function frame(scene::Scene, framenumber::Int64, cpm::CellPotts)
    
    ModelStep!(cpm)

    migrationIndex = findfirst(x->x isa MigrationPenalty, cpm.penalties)
    nodeMemory = cpm.penalties[migrationIndex].nodeMemory

    background("grey90")
    setline(1.0)
    drawCells(cpm; colorby=:type, cellsize = (5,5))
    drawHeatmap(nodeMemory; colormap=:Purples, cellsize=(5,5))
    
end

animate(mymovie,
        [
            Scene(mymovie, (scene,framenumber) -> frame(scene,framenumber,cpm), 1:300)
        ],
    creategif=true,
    pathname="migration2.gif")