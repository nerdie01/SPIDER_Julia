module MHDSpectral

include("../../../functions/VizUtils.jl")
using MeshGrid, GLMakie, Random, FFTW, LinearAlgebra
using .VizUtils

include("SpectralTerms.jl")
using .SpectralTerms

export mhd_spectral
function mhd_spectral(fω::Function, fj::Function, N::Int, dt::AbstractFloat, butcher_A::AbstractMatrix, butcher_b::AbstractVector, timesteps::Int, ν::AbstractFloat, η::AbstractFloat, explicit::Bool=false, vis::Bool=true, vis_path_prefix="", forcing::Function=(x,y)->0*x)
    function rk_step_explicit!(f, ω::AbstractMatrix, A::AbstractMatrix)
        # translated from octave using an LLM
        s = length(butcher_b)
        k_ω = Vector{Matrix{ComplexF64}}(undef, s)
        k_A = Vector{Matrix{ComplexF64}}(undef, s)

        k_ω[1], k_A[1] = f(ω, A)

        for i in 2:s
            ω_stage = ω .+ dt .* sum(butcher_A[i,j] .* k_ω[j] for j in 1:i-1; init=zeros(ComplexF64, size(ω)))
            A_stage = A .+ dt .* sum(butcher_A[i,j] .* k_A[j] for j in 1:i-1; init=zeros(ComplexF64, size(A)))
            k_ω[i], k_A[i] = f(ω_stage, A_stage)
        end

        ω .+= dt .* sum(butcher_b[i] .* k_ω[i] for i in 1:s; init=zeros(ComplexF64, size(ω)))
        A .+= dt .* sum(butcher_b[i] .* k_A[i] for i in 1:s; init=zeros(ComplexF64, size(A)))
    end

    function rk_step_implicit!(f, ω::AbstractMatrix, A::AbstractMatrix, max_iter::Int=1000, tol::Float64=1e-6)
        # also LLM translated. I think it corrected some errors too, cool!
        s = length(butcher_b)

        # Initial guess: evaluate RHS at current (ω, A) and reuse for all stages
        k0_ω, k0_A = f(ω, A)
        k_ω = [copy(k0_ω) for _ in 1:s]
        k_A = [copy(k0_A) for _ in 1:s]

        converged = false
        for iter in 1:max_iter
            # Compute stage values from current K estimates
            ω_stages = [ω .+ dt .* sum(butcher_A[j,l] .* k_ω[l]
                            for l in 1:s; init=zeros(ComplexF64, size(ω)))
                        for j in 1:s]
            A_stages = [A .+ dt .* sum(butcher_A[j,l] .* k_A[l]
                            for l in 1:s; init=zeros(ComplexF64, size(A)))
                        for j in 1:s]

            # Re-evaluate RHS at each stage
            k_ω_new = Vector{Matrix{ComplexF64}}(undef, s)
            k_A_new = Vector{Matrix{ComplexF64}}(undef, s)
            for j in 1:s
                k_ω_new[j], k_A_new[j] = f(ω_stages[j], A_stages[j])
            end

            # Relative Frobenius norm residual, averaged over stages and fields
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

        # Accumulate weighted stage contributions in-place
        ω .+= dt .* sum(butcher_b[j] .* k_ω[j]
                        for j in 1:s; init=zeros(ComplexF64, size(ω)))
        A .+= dt .* sum(butcher_b[j] .* k_A[j]
                        for j in 1:s; init=zeros(ComplexF64, size(A)))
    end

    x, y = meshgrid(LinRange(0, 2π, N+1)[1:N], LinRange(0, 2π, N+1)[1:N])

    u = zeros((N,N,timesteps))
    v = zeros((N,N,timesteps))
    Bx = zeros((N,N,timesteps))
    By = zeros((N,N,timesteps))
    p = zeros((N,N,timesteps))

    kvec = collect(0:N-1)
    kvec[kvec .> N÷2] .-= N

    kx = reshape(kvec, 1, N)
    ky = reshape(kvec, N, 1)

    mask1d = abs.(kvec) .<= N÷3
    mask = mask1d' .* mask1d

    ∇2 = -kx.^2 .- ky.^2
    ∇2[1,1] = 1.0

    ∂x = im.*kx .* mask
    ∂y = im.*ky .* mask
    inv_∇2 = complex((1.0 ./ ∇2) .* mask)

    ∂x[1,1] = 0
    ∂y[1,1] = 0
    inv_∇2[1,1] = 0

    kx = kx .* mask
    ky = ky .* mask

    ω = fft(fω.(x,y)) .* mask
    A = -inv_∇2 .* fft(fj.(x,y)) .* mask
    
    println("=== starting spectral ===")
    for t ∈ 1:timesteps
        uhat = -∂y.*inv_∇2.*ω
        vhat = ∂x.*inv_∇2.*ω
        Bxhat = ∂y.*A
        Byhat = -∂x.*A

        u[:,:,t] = real(ifft(uhat))
        v[:,:,t] = real(ifft(vhat))
        Bx[:,:,t] = real(ifft(Bxhat))
        By[:,:,t] = real(ifft(Byhat))

        ux = real(ifft(∂x .* uhat))
        uy = real(ifft(∂y .* uhat))
        vx = real(ifft(∂x .* vhat))
        vy = real(ifft(∂y .* vhat))

        Bxx = real(ifft(∂x .* Bxhat))
        Bxy = real(ifft(∂y .* Bxhat))
        Byx = real(ifft(∂x .* Byhat))
        Byy = real(ifft(∂y .* Byhat))

        p_rhs = Bxx.^2 .+ Byy.^2 .+ 2 .*Bxy.*Byx .- ux.^2 .- vy.^2 .- 2 .*uy.*vx
        p[:,:,t] = real(ifft(inv_∇2 .* fft(p_rhs) .* mask))

        rhs = (ω_, A_) -> SpectralTerms.mhd_terms(ω_, A_, ∂x, ∂y, uhat, vhat, ν, η, forcing(x,y), mask)
        if explicit
            rk_step_explicit!(rhs, ω, A)
            continue
        end

        rk_step_implicit!(rhs, ω, A)

        println(t)
    end

    if vis
        println("=== visualizing spectral ===")
        VizUtils.quickanim(p, timesteps, 10, vis_path_prefix * "p.gif")
        VizUtils.quickanim(Bx, timesteps, 10, vis_path_prefix * "Bx.gif")
        VizUtils.quickanim(By, timesteps, 10, vis_path_prefix * "By.gif")
        VizUtils.quickanim(u, timesteps, 10, vis_path_prefix * "u.gif")
        VizUtils.quickanim(v, timesteps, 10, vis_path_prefix * "v.gif")
    end

    return (x[1,:],x[1,:],u,v,Bx,By,p)
end

end