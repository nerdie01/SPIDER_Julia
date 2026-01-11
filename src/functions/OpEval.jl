module OpEval

op_funcs::Dict{String, Function} = Dict([
    ("" => x -> x),
    ("^2" => x -> x^2),
    (".^2" => x -> x.^2),
    ("sin" => x -> sin(x)),
    ("sin." => x -> sin.(x)),
    ("cos" => x -> cos(x)),
    ("cos." => x -> cos.(x)),
    ("tan" => x -> tan(x)),
    ("tan." => x -> tan.(x)),
    ("exp" => x -> exp(x)),
    ("exp." => x -> exp.(x)),
    ("+" => (x,y) -> x + y),
    (".+" => (x,y) -> x .+ y),
    ("-" => (x,y) -> x - y),
    (".-" => (x,y) -> x .- y),
    ("*" => (x,y) -> x * y),
    (".*" => (x,y) -> x .* y),
    ("^n" => (x, y) -> x^y),
    (".^n" => (x, y) -> x.^y)
])

export op_eval
# this could be improved, but idk how to dynamically pass variables into functions... oh well works for now
function op_eval(field::AbstractDict, ops::AbstractDict)
    result = Nothing
    mod_field = copy(field)
    for (key::String, value::AbstractArray{String}) in ops
        if result != Nothing
            mod_field["res"] = result
        end
        if only(methods(op_funcs[key])).nargs - 1 == 1
            result = op_funcs[key](mod_field[value[1]])
        else
            result = op_funcs[key](mod_field[value[1]], mod_field[value[2]])
        end
    end

    return result
end

end

