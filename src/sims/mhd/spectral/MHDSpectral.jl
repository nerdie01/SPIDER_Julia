module MHDSpectral

include("../../../functions/VizUtils.jl")
using MeshGrid, GLMakie, Random, FFTW
using .VizUtils

include("../../CommonSimTools.jl")
using .CommonSimTools

include("SpectralTerms.jl")
using .SpectralTerms

export mhd_spectral
function mhd_spectral(fω::Function, fj::Function, N::Int, dt::AbstractFloat, butcher_A::AbstractMatrix, butcher_b::AbstractVector, timesteps::Int, ν::AbstractFloat, η::AbstractFloat, vis::Bool=true, vis_path_prefix="", forcing::Function=(x,y)->0*x, explicit::Bool=false)
    function rk_step_explicit(f, ω::Matrix{ComplexF64}, A::Matrix{ComplexF64}, dt::Float64, butcher::ButcherTableau)::Tuple{Matrix{ComplexF64}, Matrix{ComplexF64}}
        # translated from octave using an LLM
        s = length(butcher.B)
        k_ω = Vector{Matrix{ComplexF64}}(undef, s)
        k_A = Vector{Matrix{ComplexF64}}(undef, s)

        k_ω[1], k_A[1] = f(ω, A)

        for i in 2:s
            ω_stage = ω .+ dt .* sum(butcher.A[i,j] .* k_ω[j] for j in 1:i-1; init=zeros(ComplexF64, size(ω)))
            A_stage = A .+ dt .* sum(butcher.A[i,j] .* k_A[j] for j in 1:i-1; init=zeros(ComplexF64, size(A)))
            k_ω[i], k_A[i] = f(ω_stage, A_stage)
        end

        ω_new = ω .+ dt .* sum(butcher.B[i] .* k_ω[i] for i in 1:s; init=zeros(ComplexF64, size(ω)))
        A_new = A .+ dt .* sum(butcher.B[i] .* k_A[i] for i in 1:s; init=zeros(ComplexF64, size(A)))

        return ω_new, A_new
    end

    butcher::ButcherTableau = ButcherTableau(butcher_A, butcher_b)

    x::Matrix{Float64}, y::Matrix{Float64} = meshgrid(LinRange(0, 2π, N+1)[1:N], LinRange(0, 2π, N+1)[1:N])

    u::Array{Float64} = zeros((N,N,timesteps))
    v::Array{Float64} = zeros((N,N,timesteps))
    Bx::Array{Float64} = zeros((N,N,timesteps))
    By::Array{Float64} = zeros((N,N,timesteps))
    p::Array{Float64} = zeros((N,N,timesteps))

    kvec::Vector{Int} = collect(0:N-1)
    kvec[kvec .> N÷2] .-= N

    kx::Array{ComplexF64} = reshape(kvec, 1, N)
    ky::Array{ComplexF64} = reshape(kvec, N, 1)

    mask1d::Vector{ComplexF64} = abs.(kvec) .<= N÷3
    mask::Array{ComplexF64} = mask1d' .* mask1d

    ∇2::Array{ComplexF64} = -kx.^2 .- ky.^2
    ∇2[1,1] = 1.0

    ∂x::Matrix{ComplexF64} = im.*kx .* mask
    ∂y::Matrix{ComplexF64} = im.*ky .* mask
    inv_∇2::Matrix{ComplexF64} = complex((1.0 ./ ∇2) .* mask)

    ∂x[1,1] = 0
    ∂y[1,1] = 0
    inv_∇2[1,1] = 0

    kx = kx .* mask
    ky = ky .* mask

    ω::Array{ComplexF64} = fft(fω.(x,y)) .* mask
    A::Array{ComplexF64} = -inv_∇2 .* fft(fj.(x,y)) .* mask
    
    for t ∈ 1:timesteps
        uhat::Matrix{ComplexF64} = -∂y.*inv_∇2.*ω
        vhat::Matrix{ComplexF64} = ∂x.*inv_∇2.*ω
        Bxhat::Matrix{ComplexF64} = ∂y.*A
        Byhat::Matrix{ComplexF64} = -∂x.*A

        u[:,:,t] = real(ifft(uhat))
        v[:,:,t] = real(ifft(vhat))
        Bx[:,:,t] = real(ifft(Bxhat))
        By[:,:,t] = real(ifft(Byhat))

        ux::Matrix{ComplexF64} = real(ifft(∂x .* uhat))
        uy::Matrix{ComplexF64} = real(ifft(∂y .* uhat))
        vx::Matrix{ComplexF64} = real(ifft(∂x .* vhat))
        vy::Matrix{ComplexF64} = real(ifft(∂y .* vhat))

        Bxx::Matrix{ComplexF64} = real(ifft(∂x .* Bxhat))
        Bxy::Matrix{ComplexF64} = real(ifft(∂y .* Bxhat))
        Byx::Matrix{ComplexF64} = real(ifft(∂x .* Byhat))
        Byy::Matrix{ComplexF64} = real(ifft(∂y .* Byhat))

        p_rhs = Bxx.^2 .+ Byy.^2 .+ 2 .*Bxy.*Byx .- ux.^2 .- vy.^2 .- 2 .*uy.*vx
        p[:,:,t] = real(ifft(inv_∇2 .* fft(p_rhs) .* mask))

        rhs = (ω_, A_) -> SpectralTerms.mhd_terms(ω_, A_, ∂x, ∂y, uhat, vhat, ν, η, forcing(x,y), mask)

        ω, A = rk_step_explicit(rhs, ω, A, dt, butcher)
    end

    if vis
        VizUtils.quickanim(p, timesteps, 10, vis_path_prefix * "p.gif")
        VizUtils.quickanim(Bx, timesteps, 10, vis_path_prefix * "Bx.gif")
        VizUtils.quickanim(By, timesteps, 10, vis_path_prefix * "By.gif")
        VizUtils.quickanim(u, timesteps, 10, vis_path_prefix * "u.gif")
        VizUtils.quickanim(v, timesteps, 10, vis_path_prefix * "v.gif")
    end
end

end