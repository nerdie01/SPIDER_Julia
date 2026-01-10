module ReportModel

using Printf

function report_identified_model(cs::AbstractMatrix{<:Real}, residuals::AbstractVector{<:Real}, scales::AbstractVector{<:Real}, labels::AbstractVector{<:AbstractString}, gamma::Real, hide_coef_tol::Real=0.005)::Int

    ratio::Vector{Float64} = (residuals[1:end-1]) ./ (residuals[2:end])
    jumps::Vector{Int} = findall(ratio .> gamma)
    k::Int = jumps[end] + 1

    c::Vector{Float64} = cs[:, k] ./ scales
    I::Vector{Int64} = findall(x -> x ≉ 0, c)
    c = c[I]

    max_abs::Float64 = maximum(abs.(c))
    c = c ./ max_abs

    println(c)
    println(labels[I])
    println('\n',"="^20,'\n')

    if c[1] < 0
        c = -c
    end

    coef::String = ""
    for i in eachindex(I)
        if abs(abs(c[i]) - 1) < hide_coef_tol
            global coef = @sprintf("")
        else
            global coef = @sprintf("%.2g", abs(c[i]))
        end
        
        print(coef * labels[I[i]])
        if i < length(I)
            if c[i+1] > 0
                print(" + ")
            else
                print(" - ")
            end
        end
    end
    println(" = 0")

    println('\n',"="^20,'\n')
    
    return k
end

end