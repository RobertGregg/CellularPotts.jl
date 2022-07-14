using Random, Distributions
using CairoMakie
using Makie.GeometryBasics

#How large is the image?
imSize = (100,100)

#Create a random binary matrix where you can control the number ones
im = rand(Bernoulli(0.01),imSize)

#Where are those ones?
pointLoc = findall(im)

#Create another matrix to caluclate nearest distances
distances = zeros(Int,imSize)

#Distance function for two CartesianIndices
d(x::CartesianIndex, y::CartesianIndex) = sum(x->x^2, x.I .- y.I)


#Loop through all the points and save the closest distances
function calculateMinDistance!(distances,pointLoc)
    for I in CartesianIndices(distances)
        bestDistance = Inf 
        for (i,J) in enumerate(pointLoc)
            currentDistance = d(I,J)
            if i==1 || currentDistance < bestDistance
                distances[I] = i
                bestDistance = currentDistance
            end
        end
    end
end

calculateMinDistance!(distances,pointLoc)


density(distances[:], color = (:blue, 0.3),
    strokecolor = :blue, strokewidth = 3, strokearound = true)

heatmap(distances, colormap=:greys)


f = Figure()
ax = Axis(
    f[1, 1],
    xgridvisible = false,
    ygridvisible = false)

ps = [Polygon(rand(Point2f, 3) .+ Point2f(i, j))
    for i in 1:5 for j in 1:10]

poly!(ax,ps, color = :white, strokecolor = :black,strokewidth = 1, grid = false)

display(f)



#####################
using VoronoiCells
using GeometryBasics
using CairoMakie
using Makie.GeometryBasics
using Random


rng = Random.MersenneTwister(1337)
rect = Rectangle(Point2(0, 0), Point2(1, 1))
points = [Point2(rand(rng), rand(rng)) for _ in 1:1000]


tess = voronoicells(points, rect)


ps = Polygon.(tess.Cells);



f = Figure()
ax = Axis(
    f[1, 1],
    xgridvisible = false,
    ygridvisible = false)

tightlimits!(ax)
hidedecorations!(ax)

poly!(ax,ps, color = :white, strokecolor = :black,strokewidth = 1, grid = false)

display(f)