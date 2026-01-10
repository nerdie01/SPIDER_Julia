module GreedyRegression

using LinearAlgebra

function greedy_regression(G::AbstractMatrix{<:Real})::Tuple{Matrix{Float64}, Vector{Float64}, Matrix{Float64}}
    m::Int, n::Int = size(G)

    fbounds::Matrix{Float64} = zeros(n, 2)
    cs::Matrix{Float64} = zeros(n, n)
    R = Matrix(qr(G).R)
    G_rot = R[1:n, 1:n]
    
    A0::Matrix{Float64} = G_rot

    I::BitVector = trues(n)
    residuals::Vector{Float64} = zeros(n)
    k::Int = n

    while k > 0
        current_G = A0[:, I]

        F_svd = svd(current_G)
        U = F_svd.U
        S_vec = F_svd.S
        V = F_svd.V

        cs[I, k] = V[:, end]
        residuals[k] = S_vec[end]

        if k == 1
            break
        end

        candidates::Vector{Float64} = zeros(k)
        s1::Float64 = S_vec[end]
        s2::Float64 = S_vec[end-1]
        r::Float64 = s1 / s2

        for i in 1:k
            a = current_G[:, i]
            alpha::Float64 = 1.0 / norm(a)
            w = alpha * (U' * a)

            function f(Σ::Float64)
                Σ_s = Σ / s2
                term1 = (r^2 - Σ_s^2) * (s2^2 - Σ^2) * alpha^2 / (r^2 - 1)
                term2 = -w[end]^2 * (1 - Σ_s^2) / (r^2 - 1)
                term3 = -w[end-1]^2 * (r^2 - Σ_s^2) / (r^2 - 1)
                term4 = 0.0
                
                if !isempty(S_vec[1:end-2])
                    denom = (S_vec[1:end-2] ./ s2).^2 .- Σ_s^2
                    term4 = -sum((w[1:end-2].^2) ./ denom) * (r^2 - Σ_s^2) * (1 - Σ_s^2) / (r^2 - 1)
                end
                return term1 + term2 + term3 + term4
            end

            fbounds[i, 1] = f(s1)
            fbounds[i, 2] = f(s2)

            maxit::Int = 128
            threshold::Float64 = 1e-130
            g_val::Float64 = 0.0
            lb = s1
            ub = s2
            
            for _ in 1:maxit
                g_val = (lb + ub)/ 2.0
                fg = f(g_val)
                if abs(fg) < threshold
                    break
                end
                if fg > 0
                    lb = g_val
                else
                    ub = g_val
                end
            end
            candidates[i] = g_val
        end

        active_indices = findall(I)
        idx_to_remove = active_indices[argmin(candidates[1:k])]
        I[idx_to_remove] = false
        k -= 1
    end

    residuals = residuals ./ sqrt(m)
    return cs, residuals, fbounds
end

end