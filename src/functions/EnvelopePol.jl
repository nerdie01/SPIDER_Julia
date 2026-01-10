module EnvelopePol

function envelope_pol(m::Int, n::Int)::Matrix{Float64}
    deg::Int = 2*m
    coeffs_asc::Vector{Float64} = zeros(deg+1)

    for k in 0:m
        coeffs_asc[2*k+1] = binomial(m, k) * ((-1)^k)
    end

    pol::Matrix{Float64} = repeat(reshape(coeffs_asc, :, 1), 1, n)
    
    return pol
end

end