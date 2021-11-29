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