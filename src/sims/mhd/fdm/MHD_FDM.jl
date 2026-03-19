module MHD_FDM

include("../../../functions/VizUtils.jl")
using MeshGrid, GLMakie, Random, LinearAlgebra
using .VizUtils

export mhd_fdm
function mhd_fdm(fω::Function, fj::Function, N::Int, dt::AbstractFloat, butcher_A::AbstractMatrix, butcher_b::AbstractVector, timesteps::Int, ν::AbstractFloat, η::AbstractFloat, explicit::Bool=false, vis::Bool=true, vis_path_prefix="", forcing::Function=(x,y)->0*x)
    x, y::Matrix{Float64} = meshgrid(LinRange(0, 2π, N+1)[1:N], LinRange(0, 2π, N+1)[1:N])
    dx::Float64 = x[1,2] - x[1,1]

    function inv_∇2(f::AbstractMatrix{Float64}, ϕ0::AbstractMatrix{Float64}=zeros(N,N), Ω::Float64=1.0, iter::Int=1000)::Matrix{Float64}
        ϕ = copy(ϕ0)
        dx2::Float64 = dx^2

        for t in 1:iter
            for i in 1:N
                i_plus = mod1(i+1, N)
                i_minus = mod1(i-1, N)
                
                for j in 1:N
                    j_plus = mod1(j+1, N)
                    j_minus = mod1(j-1, N)

                    ϕ[i, j] = (1.0 - Ω)*ϕ[i, j] + 0.25Ω*(ϕ[i_plus, j] + ϕ[i_minus, j] + ϕ[i, j_plus] + ϕ[i, j_minus] - dx2 * f[i, j])
                end
            end
        end

        return ϕ
    end

    function ∇2(f::AbstractMatrix{Float64})::Matrix{Float64}
        plus_x = mod1.(2:N+1, N); minus_x = mod1.(0:N-1, N)
        plus_y = mod1.(2:N+1, N); minus_y = mod1.(0:N-1, N)
        
        return (f[plus_x,:] + f[minus_x,:] + f[:,plus_y] + f[:,minus_y] - 4f) ./ dx^2
    end

    function ∂x(f::AbstractMatrix{Float64})
        plus = mod1.(2:N+1, N)
        minus = mod1.(0:N-1, N)

        return (f[:,plus] - f[:,minus]) ./ 2dx
    end

    function ∂y(f::AbstractMatrix{Float64})
        plus = mod1.(2:N+1, N)
        minus = mod1.(0:N-1, N)

        return (f[plus,:] - f[minus,:]) ./ 2dx
    end

    ω = fω.(x,y)
    A = inv_∇2(-fj.(x,y))
    frc = forcing(x,y)

    u = zeros((N,N,timesteps))
    v = zeros((N,N,timesteps))
    Bx = zeros((N,N,timesteps))
    By = zeros((N,N,timesteps))
    p = zeros((N,N,timesteps))

    function fdm_terms(ω::AbstractMatrix, A::AbstractMatrix)
        ∂xω = ∂x(ω)
        ∂yω = ∂y(ω)
        u_  = inv_∇2(-∂yω)
        v_  = inv_∇2(∂xω)
        curr_  = -∇2(A)
        Bx_ =  ∂y(A)
        By_ = -∂x(A)

        adv_ω = u_ .* ∂xω .+ v_ .* ∂yω
        lorentz = Bx_ .* ∂x(curr_) .+ By_ .* ∂y(curr_)
        visc = ν .* ∇2(ω) .+ frc
        nl_ω = -adv_ω .+ lorentz .+ visc

        nl_A = -(u_ .* ∂x(A) .+ v_ .* ∂y(A)) .+ η .* ∇2(A)

        return (nl_ω, nl_A)
    end

    # it looks like this method is usually too unstable to be useful at all
    # a timestep of 2^-6 for spectral matches up with ~2^-15 for finite difference, and it doesn't seem to be stability related
    # idk what's going on but i'll just handle explicit schemes implicitly too
    function rk_step_explicit!(ω::AbstractMatrix, A::AbstractMatrix)
        s = length(butcher_b)
        k_ω = Vector{Matrix{Float64}}(undef, s)
        k_A = Vector{Matrix{Float64}}(undef, s)

        k_ω[1], k_A[1] = fdm_terms(ω, A)
        for i ∈ 2:s
            ω_s = ω .+ dt .* sum(butcher_A[i,j] .* k_ω[j] for j ∈ 1:i-1; init=zeros(N, N))
            A_s = A .+ dt .* sum(butcher_A[i,j] .* k_A[j] for j ∈ 1:i-1; init=zeros(N,N))
            k_ω[i], k_A[i] = fdm_terms(ω_s, A_s)
        end

        ω .+= dt .* sum(butcher_b[i] .* k_ω[i] for i ∈ 1:s; init=zeros(N,N))
        A .+= dt .* sum(butcher_b[i] .* k_A[i] for i ∈ 1:s; init=zeros(N,N))
    end

    function rk_step_implicit!(ω::AbstractMatrix, A::AbstractMatrix, max_iter::Int=1000, tol::AbstractFloat=1e-6)
        s = length(butcher_b)

        k0_ω, k0_A = fdm_terms(ω, A)
        k_ω = [copy(k0_ω) for _ in 1:s]
        k_A = [copy(k0_A) for _ in 1:s]

        converged = false
        for _ ∈ 1:max_iter
            ω_stages = [ω .+ dt .* sum(butcher_A[j,l] .* k_ω[l] for l ∈ 1:s; init=zeros(N, N)) for j ∈ 1:s]
            A_stages = [A .+ dt .* sum(butcher_A[j,l] .* k_A[l] for l ∈ 1:s; init=zeros(N, N)) for j ∈ 1:s]

            k_ω_new = Vector{Matrix{Float64}}(undef, s)
            k_A_new = Vector{Matrix{Float64}}(undef, s)

            for j ∈ 1:s
                k_ω_new[j], k_A_new[j] = fdm_terms(ω_stages[j], A_stages[j])
            end

            # breaks when sufficiently close
            err = sum(1:s) do j
                norm(k_ω_new[j] .- k_ω[j]) / max(norm(k_ω[j]), 1e-12) +
                norm(k_A_new[j] .- k_A[j]) / max(norm(k_A[j]), 1e-12)
            end / (2s)

            k_ω = k_ω_new
            k_A = k_A_new

            if err < tol
                converged = true
                break
            end
        end

        converged || @warn "Fixed-point iteration did not converge" max_iter tol
        
        ω .+= dt .* sum(butcher_b[j] .* k_ω[j] for j in 1:s; init=zeros(N, N))
        A .+= dt .* sum(butcher_b[j] .* k_A[j] for j in 1:s; init=zeros(N, N))
    end

    function update_fields!(t::Int, ω::AbstractMatrix, A::AbstractMatrix)
        u[:,:,t] = inv_∇2(-∂y(ω))
        v[:,:,t] = inv_∇2(∂x(ω))
        curr = -∇2(A)
        Bx[:,:,t] = ∂y(A)
        By[:,:,t] = -∂x(A)

        ux = ∂x(u[:,:,t]); uy = ∂y(u[:,:,t])
        vx = ∂x(v[:,:,t]); vy = ∂y(v[:,:,t])
        p[:,:,t] = inv_∇2(2 .* ux .* vy .- 2 .* uy .* vx .- ∂x(curr .* By[:,:,t]) .+ ∂y(curr .* Bx[:,:,t]))
    end

    update_fields!(1, ω, A)

    println("=== starting FDM ===")
    for t ∈ 2:timesteps
        rk_step_implicit!(ω, A)
        update_fields!(t, ω, A)
        println(t)
    end

    if vis
        println("=== visualizing FDM ===")
        VizUtils.quickanim(p, timesteps, 10, vis_path_prefix * "p.gif")
        VizUtils.quickanim(Bx, timesteps, 10, vis_path_prefix * "Bx.gif")
        VizUtils.quickanim(By, timesteps, 10, vis_path_prefix * "By.gif")
        VizUtils.quickanim(u, timesteps, 10, vis_path_prefix * "u.gif")
        VizUtils.quickanim(v, timesteps, 10, vis_path_prefix * "v.gif")
    end

    return (x[1,:],x[1,:],u,v,Bx,By,p)
end

end