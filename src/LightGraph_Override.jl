using LightGraphs
import LightGraphs: induced_subgraph, SimpleGraph

using TimerOutputs

using CellularPotts
M = ModelParameters()
CPM = CellPotts(M)

g = CPM.graph.network
cellIdx = findall(isequal(1), CPM.graph.σ)


function induced_subgraph(g::T, vlist::AbstractVector{U}) where T <: AbstractGraph where U <: Integer
    allunique(vlist) || throw(ArgumentError("Vertices in subgraph list must be unique"))
    h = T(length(vlist))
    newvid = Dict{U,U}()
    vmap = Vector{U}(undef, length(vlist))
    for (i, v) in enumerate(vlist)
        newvid[v] = U(i)
        vmap[i] = v
    end

    vset = Set(vlist)
    for s in vlist
        for d in outneighbors(g, s)
            # println("s = $s, d = $d")
            if d in vset && has_edge(g, s, d)
                newe = Edge(newvid[s], newvid[d])
                add_edge!(h, newe)
            end
        end
    end
    return h, vmap
end



function testsub(g::T, vlist::AbstractVector{U}) where T <: AbstractGraph where U <: Integer
    reset_timer!()

    @timeit "section 1" begin
        allunique(vlist) || throw(ArgumentError("Vertices in subgraph list must be unique"))
        h = T(length(vlist))
        newvid = Dict{U,U}()
        vmap = Vector{U}(undef, length(vlist))
        for (i, v) in enumerate(vlist)
            newvid[v] = U(i)
            vmap[i] = v
        end 
    end

    @timeit "section 2" begin
        vset = Set(vlist)
        for s in vlist
            for d in outneighbors(g, s)
                # println("s = $s, d = $d")
                if d in vset && has_edge(g, s, d) && !has_edge(h, newvid[s], newvid[d])
                    newe = Edge(newvid[s], newvid[d])
                    add_edge!(h, newe)
                end
            end
        end
    end
    print_timer()
    return h
end



function testsub(g::SimpleGraph, cellIdx::AbstractVector{Int})
    
    fadjlist = copy.(g.fadjlist[cellIdx])
    #vDict = Dict(cellIdx .=> eachindex(cellIdx))

    for i in eachindex(fadjlist)
        filter!(x -> x ∈ cellIdx, fadjlist[i])

        for j in eachindex(fadjlist[i])
            fadjlist[i][j] = searchsortedfirst(cellIdx,fadjlist[i][j])
        end
    end

   return SimpleGraph( sum(length,fadjlist) ÷ 2, fadjlist)
end


function testsub_timed(g::SimpleGraph, cellIdx::AbstractVector{Int})

    reset_timer!()
    
    @timeit "make fadj list" fadjlist = copy.(g.fadjlist[cellIdx])
    #@timeit "create dict" vDict = Dict(cellIdx .=> eachindex(cellIdx))

    for i in eachindex(fadjlist)
        @timeit "filter" filter!(x -> x ∈ cellIdx, fadjlist[i])

        for j in eachindex(fadjlist[i])
            #@timeit "map dict" fadjlist[i][j] = vDict[ fadjlist[i][j] ]
            @timeit "findfirst" fadjlist[i][j] = findfirst(isequal(fadjlist[i][j]), cellIdx)
        end
    end

    @timeit "create graph" gsub =  SimpleGraph( sum(length,fadjlist) ÷ 2, fadjlist)

    print_timer()
    return gsub
end


using BenchmarkTools

@benchmark g[cellIdx]

@benchmark testsub($g,$cellIdx)


@benchmark vDict = Dict(cellIdx .=> eachindex(cellIdx))

#Creating a view of a graph?

g = SimpleGraph(CellularPotts.genAdj(4,4, true))
cellIdx = [1,4,5,13,16]


function testsub(g::SimpleGraph, cellIdx::AbstractVector{Int})

    gSubRaw = view(g.fadjlist, cellIdx)
    fadjlist = filter.(x -> x ∈ cellIdx, gSubRaw)

    for list in fadjlist
        for (j,vertex) in enumerate(list)
            list[j] = searchsortedfirst(cellIdx,vertex)
        end
    end

    return SimpleGraph( sum(length,fadjlist) ÷ 2, fadjlist)
end



function testsub(g::SimpleGraph, cellIdx::AbstractVector{Int})

    gSubRaw = view(g.fadjlist, cellIdx)
    fadjlist = filter.(x -> x ∈ cellIdx, gSubRaw)
    d = Dict(cellIdx .=> eachindex(cellIdx))

    for list in fadjlist
        for (j,vertex) in enumerate(list)
            #list[j] = searchsortedfirst(cellIdx,vertex)
            list[j] = d[vertex]
        end
    end

    return SimpleGraph( sum(length,fadjlist) ÷ 2, fadjlist)
end
