#Using discrete amounts for molecules might be better than ODEs?

using Catalyst
sir_model = @reaction_network begin
    β, S + I --> 2I
    ν, I --> R
end β ν


p     = (0.1/1000,0.01)   
u₀    = [999,1,0]
tspan = (0.0,250.0)
prob  = DiscreteProblem(sir_model, u₀, tspan, p)


using DiffEqJumps
jump_prob = JumpProblem(sir_model, prob, Direct())


sol = solve(jump_prob, SSAStepper())


using Plots
plot(sol)