
abstract type AbstractCell end

struct Cell <: AbstractCell
    name::Symbol
    id::Int
    volume::Int
end

macro newcell(name, base, fields)
    base_type = Core.eval(@__MODULE__, base)
    base_fieldnames = fieldnames(base_type)
    base_types = [t for t in base_type.types]
    base_fields = [:($f::$T) for (f, T) in zip(base_fieldnames, base_types)]
    res = :(mutable struct $(esc(name)) <: AbstractCell end)
    push!(res.args[end].args, base_fields...)
    push!(res.args[end].args, map(esc, fields.args)...)
    return res
end


@newcell Macrophage Cell begin
    IL6::Float64
    isHappy::Bool
end

m = Macrophage(:name, 1, 10, 52.31, true)

t = Macrophage