module PickSubdomains

using Random

function pick_subdomains(size_of_data::Vector{Int}, size_vec::AbstractVector{Int}, buffer::Int, nw::Int, seed::Int)::Matrix{Int}
    Random.seed!(seed)

    dim::Int = length(size_of_data)

    corners::Matrix{Int} = zeros(Int, dim, nw)
    for d in 1:dim
        max_val::Int = size_of_data[d] - size_vec[d] - 2*buffer
        corners[d, :] = rand(1:max_val, nw) .+ buffer
    end
    return corners
end

end