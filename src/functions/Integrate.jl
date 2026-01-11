module Integrate

using LinearAlgebra

export poly_integrate
function poly_integrate(data::AbstractArray, derivs::AbstractVector{Int}, grid::AbstractVector{<:AbstractVector}, corners::AbstractMatrix{Int}, size_vec::AbstractVector{Int}, pol::AbstractMatrix)

    n = size(corners, 1)
    num_windows = size(corners, 2)

    data0 = dropdims(data, dims=Tuple(findall(==(1), size(data))))

    deriv_counts = zeros(Int, n)
    for d in derivs
        deriv_counts[d] += 1
    end

    pol0 = copy(pol)
    poly_sign = 1
    derivs0 = copy(derivs)
    l = size(pol0, 1)

    while !isempty(derivs0)
        d = popfirst!(derivs0)
        poly_sign = -poly_sign
        for ll in 2:l
            pol0[ll-1, d] = (ll-1) * pol0[ll, d]
        end
        pol0[l, d] = 0.0
    end

    vals = zeros(eltype(data0), num_windows)

    for i in 1:num_windows
        specific_subsample = Vector{UnitRange{Int}}(undef, n)
        for j in 1:n
            specific_subsample[j] = corners[j, i]:(corners[j,i] + size_vec[j] - 1)
        end

        conversion_factor = poly_sign
        small_data = data0[specific_subsample...]

        for j in 1:n
            x_spec = grid[j][specific_subsample[j]]

            m = 2 / (x_spec[end] - x_spec[1])
            b = -(x_spec[end] + x_spec[1]) / (x_spec[end] - x_spec[1])
            y = m .* x_spec .+ b

            conversion_factor *= m^deriv_counts[j]
            poly_weight = [evalpoly(y_i, pol0[:, j]) for y_i in y]

            # trapezoidal rule weights
            N = length(x_spec)
            hl = zeros(N)
            hr = zeros(N)
            if N > 1
                hl[2:end] = @views x_spec[2:end] - x_spec[1:end-1]
                hr[1:end-1] = @views x_spec[2:end] - x_spec[1:end-1]
            end

            weights = poly_weight .* (hl .+ hr) ./ 2
            weights_shape = (N, ntuple(_ -> 1, ndims(small_data)-1)...)
            weighted_data = small_data .* reshape(weights, weights_shape)
            small_data = dropdims(sum(weighted_data, dims=1), dims=1)
        end

        vals[i] = small_data[] * conversion_factor
    end

    return vals
end

end