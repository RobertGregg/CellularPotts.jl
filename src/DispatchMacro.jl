#Want to replicate this
function f(m::model, ss::NamedTuple{(:s1, :s2, :s3, :s4...), Tuple{S1, S2, S3, S4...}})

    T = fieldtypes(typeof(ss))
    total = 0
    for s in values(ss)
        if s isa T[1]
            total += f(m,s)
        elseif s isa T[2]
            total += f(m,s)
        elseif s isa T[3]
            total += f(m,s)
        elseif s isa T[4]
            total += f(m,s)
            ...
        end
    end

    return total
end



function _unionsplit(types, call)
    MacroTools.@capture(call, f_(arg_, args__))
    thetypes = types.args
    first_type, rest_types = Iterators.peel(thetypes)
    code = :(if $arg isa $first_type
               $call
             end)
    the_args = code.args
    for next_type in rest_types
        clause = :(if $arg isa $next_type # use `if` so this parses, then change to `elseif`
                     $call
                   end)
        clause.head = :elseif
        push!(the_args, clause)
        the_args = clause.args
    end
    push!(the_args, call) # The last one uses Julia's dispatch system.
    return code
end