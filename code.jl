using Images
using ColorSchemes
using Random
using LoopVectorization
using BenchmarkTools

# For getting the image
function imagesc(
    data::AbstractMatrix{<:Real};
    colorscheme::ColorScheme=ColorSchemes.viridis,
    maxsize::Integer=512, rangescale=:extrema
    )

    s = maximum(size(data))
    if s > maxsize
    return imagesc(imresize(data, ratio=maxsize/s);   # imresize from Images.jl
            colorscheme, maxsize, rangescale)
    end
    return get(colorscheme, data, rangescale) # get(...) from ColorSchemes.jl
end

# Coordinates to index and vice versa
function idxToCoord(idx, size)
    return (mod(idx-1, size) + 1, div(idx-1, size) + 1)
end
function coordToIdx(i, j, size)
    return (j-1)*size + i
end

# Get the adjacency list for a square lattice
function getSqAdj(size, NN)
    @inline function dist(i1,j1, i2, j2)
        Float32(sqrt((i1-i2)^2 + (j1-j2)^2))
    end

    @inline function getweight(dr)
        # Some function of dr
        return 1f0
    end

    adj = Vector{Tuple{Int32, Float32}}[Tuple{Int32,Float32}[] for i in 1:(size*size)]

    for idx in 1:(size*size)
        vert_i, vert_j = idxToCoord(idx, size)
        for dj in (-NN):NN
            for di in (-NN):NN
                # Include self connection?
                if di == 0 && dj == 0
                    continue
                end

                # Periodicity
                conn_i = vert_i + di > size ? vert_i + di - size : vert_i + di
                conn_i = conn_i < 1 ? conn_i + size : conn_i

                conn_j = vert_j + dj > size ? vert_j + dj - size : vert_j + dj
                conn_j = conn_j < 1 ? conn_j + size : conn_j

                weight = getweight(dist(vert_i, vert_j, conn_i, conn_j))
                if weight != 0
                    push!(adj[idx], (coordToIdx(conn_i, conn_j, size), weight))
                end
            end
        end
    end

    # Sort so accesses are in order
    for i in 1:length(adj)
        sort!(adj[i], by = x -> x[1])
    end

    return adj
end

# All the connection idxs and weights for a vertex
struct Connections
    idxs::Vector{Int32}
    weights::Vector{Float32}

    Connections() = new(Int32[], Float32[])
end

@inline Base.eachindex(conns::Connections) = Base.eachindex(conns.idxs) 

# A list of vertex connections for a graph
struct AdjList{C} <: AbstractVector{C}
    data::Vector{C}
end

@inline Base.size(A::AdjList) = size(A.data)
@inline Base.getindex(A::AdjList, i) = A.data[i]
@inline Base.setindex!(A::AdjList, v, i) = (A.data[i] = v)
@inline Base.length(A::AdjList) = length(A.data)
@inline Base.eachindex(A::AdjList) = Base.eachindex(A.data)

AdjList(len) = AdjList{Connections}(Connections[Connections() for i in 1:len])

struct AdjListRef{C} <: AbstractVector{C}
    data::Vector{Ref{C}}
end

@inline Base.size(A::AdjListRef) = size(A.data)
@inline Base.getindex(A::AdjListRef, i) = A.data[i][]
@inline Base.setindex!(A::AdjListRef, v, i) = (A.data[i][] = v)
@inline Base.length(A::AdjListRef) = length(A.data)
@inline Base.eachindex(A::AdjListRef) = Base.eachindex(A.data)

AdjListRef(len) = AdjListRef{Connections}(Ref{Connections}[Ref{Connections}(Connections()) for i in 1:len])

function AdjListToRef(adjlist)
    adjlistref = AdjListRef(length(adjlist))
    for i in eachindex(adjlist)
        adjlistref[i] = deepcopy(adjlist[i])
    end
    return adjlistref
end

# Convert the adjacency list using tuples to the adjacency list using Connections
function adjTupToAdjList(adjtup)
    adjlist = AdjList(length(adjtup))
    for vert_idx in eachindex(adjtup)
        for tuple in adjtup[vert_idx]
            push!(adjlist[vert_idx].idxs, tuple[1])
            push!(adjlist[vert_idx].weights, tuple[2])
        end
    end
    return adjlist
end

struct Graph
    state::Vector{Float32}
    adjtup::Vector{Vector{Tuple{Int32,Float32}}}
    adjlist::AdjList{Connections}
    adjref::AdjListRef{Connections}
end

randomState(size) = 2f0 .* (rand(Float32, size*size) .- .5f0)

function Graph(size, NN = 1)
    state = randomState(size)
    adjtup = getSqAdj(size, NN)
    adjlist = adjTupToAdjList(adjtup)
    adjref = AdjListToRef(adjlist)
    return Graph(state, adjtup, adjlist, adjref)
end

mutable struct Sim
    const g::Graph
    const size::Int32
    shouldrun::Bool
    isrunning::Bool
    updates::Int64
    temp::Float32
end

function Sim(size, NN = 1)
    g = Graph(size, NN)
    return Sim(g, size, true, false, 0, 1f0)
end

function startLoop(sim, whichfield)
    g = sim.g
    state = g.state
    adj = getfield(g, whichfield)
    iterator = UnitRange{Int32}(1, length(state))

    innerloop(sim, g, state, adj, iterator)
end

@inline function getEnergyFactor(state, connections::C) where C <: Vector{Tuple{Int32,Float32}}
    energy = 0.0f0
    @inbounds @simd for weight_idx in eachindex(connections)
        conn_idx = connections[weight_idx][1]
        weight = connections[weight_idx][2]

        energy += -weight * state[conn_idx]
    end
    return energy
end

@inline function getEnergyFactor(state, connections::C) where C <: Connections
    energy = 0.0f0
    weights = connections.weights
    idxs = connections.idxs
    @turbo for weight_idx in eachindex(connections.idxs)
    # @inbounds @simd for weight_idx in eachindex(connections.idxs)
        conn_idx = idxs[weight_idx]
        weight = weights[weight_idx]
        energy += -weight * state[conn_idx]
    end
    return energy
end

Base.@propagate_inbounds function innerloop(sim, g, state, adj::C, iterator) where C
    sim.isrunning = true
    while sim.shouldrun
        idx = rand(iterator)

        connections = adj[idx]

        efactor = getEnergyFactor(state, connections)
    
        beta = 1f0/(sim.temp)
         
        oldstate = state[idx]
    
        newstate = 2f0*(rand(Float32)- .5f0)
    
        ediff = efactor*(newstate-oldstate)
    
        if (ediff < 0f0 || rand(Float32) < exp(-beta*ediff))
            @inbounds state[idx] = newstate 
        end

        sim.updates += 1

        GC.safepoint()
        # If run using includ, might hang without yielding
        # Running from REPL works just fine, why?
        # yield()
    end
    sim.isrunning = false
end

# Generate an image from the current state
function genImage(sim)
    state = sim.g.state
    size = sim.size
    return imagesc(reshape(state, size, size))
end

# Relocate the internal vectors
function relocate!(sim)
    g = sim.g
    for idx in eachindex(g.adjtup)
        push!(g.adjtup[idx], (0,0))
        deleteat!(g.adjtup[idx], length(g.adjtup[idx]))

        push!(g.adjlist[idx].idxs,0)
        push!(g.adjlist[idx].weights, 0)
        deleteat!(g.adjlist[idx].idxs, length(g.adjlist[idx].idxs))
        deleteat!(g.adjlist[idx].weights, length(g.adjlist[idx].weights))
    end
end
function localize!(sim)
    g = sim.g
    g.adjlist .= deepcopy(g.adjlist)
    g.adjtup .= deepcopy(g.adjtup)
    return
end

# Can't use sleep because of bug
function halt(seconds)
    ti = time()
    while time() - ti < seconds
    end
end

function testrun(sim, adj_symb = :adjlist, sleeptime = 2; print = true)
    # Reset the seed and state
    Random.seed!(1234)
    sim.g.state .= 2f0 .* rand(Float32, sim.size*sim.size) .- 1f0
    
    # Start the loop and sleep for an amount of time
    sim.shouldrun = true

    Threads.@spawn startLoop(sim, adj_symb)
    halt(sleeptime)

    # Stop the loop
    sim.shouldrun = false
    while sim.isrunning
        sleep(1e-6)
    end

    # Gather and return the data
    if print
        println("Did testrun for $adj_symb: ")
        println("$(sim.updates) updates in $(sleeptime) seconds.")
        display(genImage(sim))

    end
    totalupdates = sim.updates
    sim.updates = 0
    return totalupdates
end




