module AutoSpider

include(joinpath(dirname(Base.active_project()), "ClassImport.jl"))
using HDF5, JSON, LinearAlgebra
using .ClassImport: GreedyRegression, Integrate, ReportModel, PickSubdomains, EnvelopePol, OpEval, FUNCTION_PATH, HDF5_PATH, JSON_PATH

export auto_spider
function auto_spider(data_f::AbstractString, options_f::AbstractString)::Dict
    options::Dict = JSON.parse(read(joinpath(JSON_PATH, options_f), String))

    fields::Dict = Dict()
    h5open(joinpath(HDF5_PATH, data_f), "r") do file
        for (key::String, value::String) in options["fields"]
            fields[key] = read(file, value)
        end
    end

    nw::Int64 = options["number_of_windows"]
    dof::Int64 = options["degrees_of_freedom"]
    dim::Int64 = options["dimensions"]
    env_pow::Int64 = options["envelope_power"]
    size_of_data::Vector{Int} = options["size_of_data"]
    size_vec::Vector{Int64} = options["integration_size"]
    buffer::Int64 = options["buffer"]
    seed::Int64 = options["seed"]

    pol::Matrix{Float64} = EnvelopePol.envelope_pol(env_pow, dim)
    corners::Matrix{Int} = PickSubdomains.pick_subdomains(size_of_data, size_vec, buffer, nw, seed)

    grid = [vec(fields[i]) for i in options["grid"]]

    labels::Vector{String} = String[]
    scales::Vector{Float64} = Float64[]

    G::Matrix{Float64} = zeros(Float64, dof*nw, 0)

    function add_library_term(label::AbstractString)
        term::Dict = options["library_terms"][label]
        ops::Dict{String, Vector{String}} = term["term"]
        derivs::Vector{Int} = term["derivs"]

        col::Vector{Float64} = Integrate.poly_integrate(OpEval.op_eval(fields, ops), derivs, grid, corners, size_vec, pol)
        G = hcat(copy(G), col)
        push!(labels, label)
        push!(scales, term["scale"])
    end

    for (key, _) in options["library_terms"]
        add_library_term(key)
    end

    norm_vec::Vector{Float64} = Integrate.poly_integrate(zeros(Tuple(size_of_data)) .+ 1, Int64[], grid, corners, size_vec, pol)
    G ./= norm_vec
    G ./= scales'

    gamma::Float64 = options["gamma"]

    res::Dict{Int, Dict{String, Any}} = Dict()
    for i::Int = 1:options["equations_to_discover"]
        res[i] = Dict()
        cs::Matrix{Float64}, residuals::Vector{Float64}, _ = GreedyRegression.greedy_regression(G)
        res[i]["k"], res[i]["label"], res[i]["coeff"], res[i]["equation"] = ReportModel.report_identified_model(cs, residuals, scales, labels, gamma)

        println(res[i]["equation"])

        res[i]["G"] = G
        res[i]["cs"] = cs
        res[i]["residuals"] = residuals

        col_norms::Vector{Float64} = [norm(G[:, j] .* cs[j, res[i]["k"]]) for j = 1:size(G, 2)]
        kill::Int = argmax(col_norms)
        G = G[:, setdiff(1:end, kill)]
        labels = labels[setdiff(1:end, kill)]
        scales = scales[setdiff(1:end, kill)]
    end

    return res
end

end