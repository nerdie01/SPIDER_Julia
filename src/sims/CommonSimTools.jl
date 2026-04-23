module CommonSimTools

export ButcherTableau
struct ButcherTableau
    A::Matrix{Float64}
    B::Vector{Float64}
end

end