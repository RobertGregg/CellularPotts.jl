using Revise
using CellularPotts

#Example Model

space = CellSpace(100,100)

initialCellState = newCellState(
    [:Epithelial, :TCells],
    [75, 50],
    [100, 10]);

initialCellState = addCellProperty(initialCellState, :isTumor, false, :Epithelial);

penalties = [
    AdhesionPenalty([0 20 20;
                    20 90 40;
                    20 40 90]),
    VolumePenalty([5,5]),
    PerimeterPenalty([5,5]),
    MigrationPenalty(10,1,[:TCells], space.gridSize)
    ]

# penalties = [
#     AdhesionPenalty([0 20 20;
#                     20 90 40;
#                     20 40 90]),
#     VolumePenalty([5,5]),
#     PerimeterPenalty([5,5])
#     ]


cpm = CellPotts(space, initialCellState, penalties);

positionCellsRandom!(cpm)

MHStep!(cpm)


#####################################################################
using GLMakie, Colors, ColorSchemes

function Edge2Grid(dim)
    gridIndices = LinearIndices(dim)

    x1 = reverse(reshape(gridIndices,dim),dims=1)'[:]
    x2 = circshift(x1,dim[2])

    y1 = reverse(reshape(reverse(gridIndices),dim),dims=2)[:]
    y2 = circshift(y1,dim[1])

    append!(x1,x1[1:dim[1]])
    append!(x2,x2[1:dim[1]])
    append!(y1,y1[1:dim[1]])
    append!(y2,y2[1:dim[1]])

    return [[id1,id2] for (id1,id2) in zip([x1;y1],[x2;y2])]
end



function plotCells(cpm::CellPotts)

    fig = Figure(resolution = (1200, 1200), backgroundcolor = RGBf0(0.98, 0.98, 0.98))
    axSim = fig[1, 1] = Axis(fig, title = "Simulation")

    cmap = ColorSchemes.nipy_spectral
    #distintCellTypes = countcelltypes(cpm) + 1

    heatmap!(axSim,
                cpm.visual,
                show_axis = false,
                colormap = :Greys_3) #cgrad(cmap, distintCellTypes, categorical=true, rev=true)
        tightlimits!.(axSim)
        hidedecorations!.(axSim) #removes axis numbers


    edgeConnectors = Edge2Grid(cpm.space.gridSize)
    (m,n) = cpm.space.gridSize

    #Generate all of the edge Connections by putting a point on each cell corner
    horizontal = [Point2f0(x, y) => Point2f0(x+1, y) for x in 0.5:m-0.5, y in 0.5:m+0.5]
    vertical = [Point2f0(x, y) => Point2f0(x, y+1) for y in 0.5:n-0.5, x in 0.5:n+0.5]
    points = vcat(horizontal[:],vertical[:])

    #Determine the transparency of the linesegments
    gridflip = rotl90(cpm.visual) #https://github.com/JuliaPlots/Makie.jl/issues/205

    #Cell borders are outlined in black
    black = RGBA{Float64}(0.0,0.0,0.0,1.0);
    clear = RGBA{Float64}(0.0,0.0,0.0,0.0);

    #Loop through all the grid connected and assign the correct color
    currentEdgeColors = [gridflip[edges[1]]==gridflip[edges[2]] ? clear : black for edges in edgeConnectors];

    linesegments!(
            axSim,
            points,
            color = currentEdgeColors,
            linewidth = 2
        )

    return fig
end


#####################################################################
using Revise
using CellularPotts

#Example Model

space = CellSpace(200,200)

initialCellState = newCellState(
    [:Epithelial],
    [500],
    [1]);

positions = [(100,100)]

initialCellState = addCellProperty(initialCellState, :positions, positions)

penalties = [
    AdhesionPenalty([0 20;
                    20 100]),
    VolumePenalty([50])
    ]


cpm = CellPotts(space, initialCellState, penalties);

positionCells!(cpm)

CellGUI(cpm)


function runModel(cpm)
    for i=1:10000
        MHStep!(cpm)
    end
end


#####################################################################
using BenchmarkTools


abstract type S end

mutable struct model 
    i::Int
end

mutable struct S1 <: S
    a::Int
end

mutable struct S2 <: S
    b::Int
end

mutable struct S3 <: S
    c::Int
end

mutable struct S4 <: S
    d::Int
end



function f(m::model, s::S1)
    return m.i + s.a
end

function f(m::model, s::S2)
    return m.i + s.b^2
end

function f(m::model, s::S3)
    return m.i + s.c^3
end

function f(m::model, s::S4)
    return m.i + s.d^4
end

m = model(1)

ss1 = Dict(:s1 => S1(5), :s2 => S2(5), :s3 => S3(5), :s4 => S4(5))

ss2 = (s1 = S1(5), s2 = S2(5), s3 = S3(5), :s4 => S4(5))

ss3 = [S1(5), S2(5), S3(5), S4(5)]

function g1(m::model, ss)
    total = 0
    for s in values(ss)
        if s isa S1
            total += f(m,s)
        elseif s isa S2
            total += f(m,s)
        elseif s isa S3
            total += f(m,s)
        elseif s isa S4
            total += f(m,s)
        end
    end

    return total
end


function g3(m::model, ss)
    total = 0

    for s in values(ss)
        total += @unionsplit((S1,S2,S3,S4), f(m,s))
    end

    return total
end

function g1(m::model, ss::T) where T<:AbstractVector
    total = 0
    @inbounds for s in ss
        total += f(m,s)
    end

    return total
end



function g2(m::model, ss)
    total = 0
    for i in keys(ss)
        total += f(m,ss[i])
    end

    return total
end


function g2(m::model, ss::T) where T<:AbstractVector
    total = 0
    @inbounds for i in eachindex(ss)
        total += f(m,ss[i])
    end

    return total
end



#####################################################################
using StatsBase, StatsPlots, Combinatorics


mRange = 0:20
neighborCount = 8

n = length(with_replacement_combinations(mRange,neighborCount))

m = zeros(n)
g = zeros(n)
g1 = zeros(Int,n)

for (i,v) = enumerate(with_replacement_combinations(mRange,neighborCount))
    m[i] = mean(v)
    g[i] = geomean(v)
    g1[i] = any(iszero,v) ? 0 : sum(v) ÷ neighborCount
end

density(m)
density!(g)


plot(
    m,
    g,
    seriestype = :histogram2d,
    c = :vik,
    nbins = 21,
    show_empty_bins = :true,
    aspect_ratio=:equal,
    xlims = (0,20),
    ylims = (0,20),
    xlabel = "Arithmetic Mean",
    ylabel = "Geometric Mean",
       )

plot!(0:20,0:20, legend=false)


plot(
    g1,
    g,
    seriestype = :histogram2d,
    c = :vik,
    nbins = 21,
    show_empty_bins = :true,
    aspect_ratio=:equal,
    xlims = (0,20),
    ylims = (0,20),
    xlabel = "My Mean",
    ylabel = "Geometric Mean",
       )

plot!(0:20,0:20, legend=false)

histogram(g, normalize=:probability)
histogram!(m, normalize=:probability)

#####################################################################