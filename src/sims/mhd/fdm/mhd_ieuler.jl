module MHD_FDM

include("../../../functions/VizUtils.jl")
using MeshGrid, GLMakie, Random
using .VizUtils

export mhd_fdm
function mhd_fdm(fω::Function, fj::Function, forcing::Function, N::Int, dt::AbstractFloat, timesteps::Int, ν::AbstractFloat, η::AbstractFloat, vis::Bool=true)
    function solve_laplace_periodic(f::AbstractMatrix{Float64}, ϕ0::AbstractMatrix{Float64}=zeros(N,N), Ω::Float64=0.5, iter::Int=1000)::Matrix{Float64}
        ϕ = copy(ϕ0)
        dx2::Float64 = dx^2

        for t in 1:iter
            for i in 1:N
                i_plus = mod1(i + 1, N)
                i_minus = mod1(i - 1, N)
                
                for j in 1:N
                    j_plus = mod1(j + 1, N)
                    j_minus = mod1(j - 1, N)

                    ϕ[i, j] = (1.0 - Ω)ϕ[i, j] + 0.25Ω*(ϕ[i_plus, j] + ϕ[i_minus, j] + ϕ[i, j_plus] + ϕ[i, j_minus] - dx2 * f[i, j])
                end
            end
        end

        return ϕ
    end

    function solve_implicit(ω_m::AbstractMatrix{Float64}, λ::Float64, ω_p0::AbstractMatrix{Float64}=zeros(N,N), Ω::Float64=0.5, iter::Int=1000)::Matrix{Float64}
        ω_p = copy(ω_p0)
        dx2::Float64 = dx^2

        for t in 1:iter
            for i in 1:N
                i_plus = mod1(i + 1, N)
                i_minus = mod1(i - 1, N)
                
                for j in 1:N
                    j_plus = mod1(j + 1, N)
                    j_minus = mod1(j - 1, N)

                    ω_p[i, j] = (1.0 - Ω)ω_p[i, j] + Ω*(λ*(ω_p[i_plus, j] + ω_p[i_minus, j] + ω_p[i, j_plus] + ω_p[i, j_minus]) + ω_m[i, j])/(1+4λ)
                end
            end
        end

        return ω_p
    end

    function laplacian(f::AbstractMatrix{Float64})::Matrix{Float64}
        plus_x = mod1.(2:N+1, N); minus_x = mod1.(0:N-1, N)
        plus_y = mod1.(2:N+1, N); minus_y = mod1.(0:N-1, N)
        
        return (f[plus_x,:] + f[minus_x,:] + f[:,plus_y] + f[:,minus_y] - 4f) ./ dx^2
    end

    function ∂x(f::AbstractMatrix{Float64})
        plus = mod1.(2:N+1, N)
        minus = mod1.(0:N-1, N)

        return (f[plus,:] - f[minus,:]) ./ 2dx
    end

    function ∂y(f::AbstractMatrix{Float64})
        plus = mod1.(2:N+1, N)
        minus = mod1.(0:N-1, N)

        return (f[:,plus] - f[:,minus]) ./ 2dx
    end

    x::Matrix{Float64}, y::Matrix{Float64} = meshgrid(LinRange(0, 2π, N), LinRange(0, 2π, N))
    dx::Float64 = x[1,2] - x[1,1]

    ω::Array{Float64} = zeros((N,N,timesteps))
    curr::Matrix{Float64} = zeros((N,N))

    ω[:,:,1] = fω.(x,y)
    curr[:,:] = fj.(x,y)

    u::Array{Float64} = zeros((N,N,timesteps))
    v::Array{Float64} = zeros((N,N,timesteps))

    u[:,:,1] = solve_laplace_periodic(∂y(ω[:,:,1]))
    v[:,:,1] = solve_laplace_periodic(-∂x(ω[:,:,1]))

    A::Array{Float64} = zeros((N,N,timesteps))
    Bx::Array{Float64} = zeros((N,N,timesteps))
    By::Array{Float64} = zeros((N,N,timesteps))

    A[:,:,1] = solve_laplace_periodic(curr[:,:])
    Bx[:,:,1] = ∂y(A[:,:,1])
    By[:,:,1] = -∂x(A[:,:,1])

    p::Array{Float64} = zeros((N,N,timesteps))

    # implicit euler loop
    for t ∈ 2:timesteps
        # potential formulation probably isn't the best anymore, should replace with Bx and By evolution equations
        adv_A::Matrix{Float64} = u[:,:,t-1].*∂x(A[:,:,t-1]) + v[:,:,t-1].*∂y(A[:,:,t-1])
        nl_A = A[:,:,t-1] .- dt*adv_A
        A[:,:,t] = solve_implicit(nl_A, η*dt/dx^2)
        Bx[:,:,t] = ∂y(A[:,:,t])
        By[:,:,t] = -∂x(A[:,:,t])

        curr = -laplacian(A[:,:,t])

        adv_ω::Matrix{Float64} = (u[:,:,t-1].*∂x(ω[:,:,t-1]) .+ v[:,:,t-1]^t.*∂y(ω[:,:,t-1]))
        lorentz::Matrix{Float64} = Bx[:,:,t].*∂x(curr) + By[:,:,t].*∂y(curr)
        nl_omega::Matrix{Float64} = ω[:,:,t-1] - dt*(adv_ω - lorentz - forcing.(x,y))

        ω[:,:,t] = solve_implicit(nl_omega, ν*dt/dx^2)
        u[:,:,t] = solve_laplace_periodic(∂y(ω[:,:,t]))
        v[:,:,t] = solve_laplace_periodic(-∂x(ω[:,:,t]))

        # i'll use a better expression for pressure when i get time
        pterm::Matrix{Float64} = u[:,:,t].*∂x(u[:,:,t]) + v[:,:,t].*∂y(v[:,:,t]) - Bx[:,:,t].*∂x(Bx[:,:,t]) - By[:,:,t].*∂y(By[:,:,t])
        p[:,:,t] = solve_laplace_periodic(∂x(pterm) + ∂y(pterm))
    end

    if vis
        VizUtils.quickanim(p, timesteps, 10, "p.gif")
        VizUtils.quickanim(A, timesteps, 10, "A.gif")
        VizUtils.quickanim(Bx, timesteps, 10, "Bx.gif")
        VizUtils.quickanim(By, timesteps, 10, "By.gif")
        VizUtils.quickanim(ω, timesteps, 10, "omega.gif")
        VizUtils.quickanim(u, timesteps, 10, "u.gif")
        VizUtils.quickanim(v, timesteps, 10, "v.gif")
    end
end

end