#Getting overlapping heatmaps
using GLMakie, AbstractPlotting, Colors

xs = LinRange(0, 10, 100)
ys = LinRange(0, 15, 100)
zs = [cos(x) * sin(y) for x in xs, y in ys]


heatmap(xs,ys,zs)


as = LinRange(0, 10, 100)
bs = LinRange(0, 15, 100)
cs = [exp(-((x-5)^2/2 + (y-7.5)^2/2)) for x in xs, y in ys]

c1 = RGBA(0,0,0,0.0);
c2 = RGBA(1,0,0,1.0);

cm = range(c1, stop=c2, length=100);

heatmap!(as,bs,cs,colormap = cm)
