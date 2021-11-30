#Using discrete amounts for molecules might be better than ODEs?
using Catalyst, DiffEqJump, Latexify, Plots

sir_model = @reaction_network begin
    β, S + I --> 2I
    ν, I --> R
end β ν


p     = (0.1/1000,0.01)   
u₀    = [999,1,0]
tspan = (0.0,250.0)
prob  = DiscreteProblem(sir_model, u₀, tspan, p)

jump_prob = JumpProblem(sir_model, prob, Direct())

sol = solve(jump_prob, SSAStepper())

plot(sol)


################################################
# Lotka-Volterra Equations
################################################

LV_model = @reaction_network begin
    α,     x --> 2x
    β̄,     x + y --> y
    δ,     x + y --> 2y
    γ,     y --> ∅
end α β̄ δ γ


# odesys = convert(ODESystem, LV_model)
# latexify(LV_model)
# latexify(odesys)

p     = (2.0, 0.02, 0.02, 1.06) #β̄ = β-δ = 0.04-0.02
u₀    = [100, 100]
tspan = (0.0,20.0)

prob  = DiscreteProblem(LV_model, u₀, tspan, p)

jump_prob = JumpProblem(LV_model, prob, RSSA())

sol = solve(jump_prob, SSAStepper())

plot(sol)
plot(sol, vars=(1,2))

######################################################
# Lotka-Volterra Equations (1-D movement on lattice)
######################################################

# Number of lattice points
N = 64

#Parameters are time and rate constant, Variables are amount of predator/prey at each point
@parameters t, h, α, β̄, δ, γ
@variables X[1:N](t) Y[1:N](t)

#Rate to jump between lattice points
rate = 1/h

rxns = [] 

for i in 1:N

    if i ≠ 1
        push!(rxns, Reaction(rate, [X[i]], [X[i-1]]))
        push!(rxns, Reaction(rate, [Y[i]], [Y[i-1]]))
    end

    if i ≠ N
        push!(rxns, Reaction(rate, [X[i]], [X[i+1]]))
        push!(rxns, Reaction(rate, [Y[i]], [Y[i+1]]))
    end

    push!(rxns, Reaction(α, [X[i]], [X[i]], [1], [2]))
    push!(rxns, Reaction(β̄, [X[i],Y[i]], [Y[i]]))
    push!(rxns, Reaction(δ, [X[i],Y[i]], [Y[i]], [1,1], [2]))
    push!(rxns, Reaction(γ, [Y[i]], nothing))

end


@named LV_spatial = ReactionSystem(rxns, t, [X;Y], [h, α, β̄, δ, γ])



p     = (1/N, 2.0, 0.02, 0.02, 1.06) #β̄ = β-δ = 0.04-0.02
u₀    = zeros(Int64, 2*N)
tspan = (0.0,20.0)

u₀[1] = 100 #100 prey in first lattice
u₀[end] = 100 #100 predators in last lattice


prob  = DiscreteProblem(LV_spatial, u₀, tspan, p)
jump_prob = JumpProblem(LV_spatial, prob, RSSA(),save_positions=(false,false))

sol = solve(jump_prob, SSAStepper(), saveat=0.01)

plot(sol, legend=nothing)


anim = @animate for currState ∈ sol.u
    bar(currState[1:N], alpha=0.5, labels="prey")
    bar!(currState[N+1:end],alpha=0.5, labels="predator")
    ylims!(0,400)
end

gif(anim, "anim.gif", fps = 60)