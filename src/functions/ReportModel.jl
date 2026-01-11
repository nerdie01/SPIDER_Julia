module ReportModel

using Printf

export report_identified_model
function report_identified_model(cs::AbstractMatrix{<:Real}, residuals::AbstractVector{<:Real}, scales::AbstractVector{<:Real}, labels::AbstractVector{<:AbstractString}, gamma::Real, hide_coef_tol::Real=0.005)::Tuple{Int, Vector{String}, Vector{Float64}, String}
    ratio::Vector{Float64} = (residuals[1:end-1]) ./ (residuals[2:end])
    jumps::Vector{Int} = findall(ratio .> gamma)
    k::Int = jumps[end] + 1

    c::Vector{Float64} = cs[:, k] ./ scales
    I::Vector{Int64} = findall(x -> x ≉ 0, c)
    c = c[I]

    max_abs::Float64 = maximum(abs.(c))
    c = c ./ max_abs

    if c[1] < 0
        c = -c
    end

    eq::String = ""
    for i in eachindex(I)
        if abs(abs(c[i]) - 1) < hide_coef_tol
            eq *= ""
        else
            eq *= repr(round(abs(c[i]); digits=3))
        end
        
        eq *= labels[I[i]]
        if i < length(I)
            if c[i+1] > 0
                eq *= " + "
            else
                eq *= " - "
            end
        end
    end
    eq *= " = 0"
    
    return (k, labels[I], c, eq)
end

end