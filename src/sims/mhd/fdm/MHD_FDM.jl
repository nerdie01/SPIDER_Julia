module MHD_FDM

include("../../../functions/VizUtils.jl")
using MeshGrid, GLMakie, Random
using .VizUtils

include("../../CommonSimTools.jl")
using .CommonSimTools

export mhd_fdm
function mhd_fdm(fω::Function, fj::Function, N::Int, dt::AbstractFloat, butcher_A::AbstractMatrix, butcher_b::AbstractVector, timesteps::Int, ν::AbstractFloat, η::AbstractFloat, vis::Bool=true, vis_path_prefix="", forcing::Function=(x,y)->0*x, explicit::Bool=false)
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

    function solve_implicit(ω_m::AbstractMatrix{Float64}, λ::Float64, ω_p0::AbstractMatrix{Float64}=zeros(N,N), Ω::Float64=1.0, iter::Int=1000)::Matrix{Float64}
        ω_p = copy(ω_p0)

        for t in 1:iter
            for i in 1:N
                i_plus = mod1(i+1, N)
                i_minus = mod1(i-1, N)
                
                for j in 1:N
                    j_plus = mod1(j+1, N)
                    j_minus = mod1(j-1, N)

                    ω_p[i, j] = (1.0 - Ω)*ω_p[i, j] + Ω*(λ*(ω_p[i_plus, j] + ω_p[i_minus, j] + ω_p[i, j_plus] + ω_p[i, j_minus]) + ω_m[i, j])/(1+4λ)
                end
            end
        end

        return ω_p
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

    butcher::ButcherTableau = ButcherTableau(butcher_A, butcher_b)

    x::Matrix{Float64}, y::Matrix{Float64} = meshgrid(LinRange(0, 2π, N+1)[1:N], LinRange(0, 2π, N+1)[1:N])
    dx::Float64 = x[1,2] - x[1,1]

    ω::Array{Float64} = fω.(x,y)

    u::Array{Float64} = zeros((N,N,timesteps))
    v::Array{Float64} = zeros((N,N,timesteps))

    u[:,:,1] = inv_∇2(-∂y(ω))
    v[:,:,1] = inv_∇2(∂x(ω))

    curr::Array{Float64} = zeros((N,N,timesteps))
    curr[:,:,1] = fj.(x,y)

    A::Array{Float64} = inv_∇2(-curr[:,:,1])

    Bx::Array{Float64} = zeros((N,N,timesteps))
    By::Array{Float64} = zeros((N,N,timesteps))
    
    Bx[:,:,1] = ∂y(A)
    By[:,:,1] = -∂x(A)

    p::Array{Float64} = zeros((N,N,timesteps))

    for t ∈ 2:timesteps
        adv_j::Matrix{Float64} = u[:,:,t-1].*∂x(curr[:,:,t-1]) .+ v[:,:,t-1].*∂y(curr[:,:,t-1])
        nl_j::Matrix{Float64}  = curr[:,:,t-1] .- dt.*adv_j
        curr[:,:,t] = solve_implicit(nl_j, η*dt/dx^2)

        A = inv_∇2(-curr[:,:,t])
        Bx[:,:,t] = ∂y(A)
        By[:,:,t] = -∂x(A)

        lorentz::Matrix{Float64} = Bx[:,:,t].*∂x(curr[:,:,t]) .+ By[:,:,t].*∂y(curr[:,:,t])
        adv_ω::Matrix{Float64}   = u[:,:,t-1].*∂x(ω) .+ v[:,:,t-1].*∂y(ω)
        nl_omega::Matrix{Float64} = ω .- dt.*(adv_ω .- lorentz .- forcing.(x,y))

        ω = solve_implicit(nl_omega, ν*dt/dx^2)

        u[:,:,t] = inv_∇2(-∂y(ω))
        v[:,:,t] = inv_∇2(∂x(ω))

        ux = ∂x(u[:,:,t]); uy = ∂y(u[:,:,t])
        vx = ∂x(v[:,:,t]); vy = ∂y(v[:,:,t])
        rhs_inertial::Matrix{Float64} = -(ux.^2 .+ 2 .*uy.*vx .+ vy.^2)
        rhs_lorentz::Matrix{Float64}  = -∂x(curr[:,:,t].*By[:,:,t]) .+ ∂y(curr[:,:,t].*Bx[:,:,t])
        p[:,:,t] = inv_∇2(rhs_inertial .+ rhs_lorentz)
    end

    if vis
        VizUtils.quickanim(p, timesteps, 10, vis_path_prefix * "p.gif")
        VizUtils.quickanim(Bx, timesteps, 10, vis_path_prefix * "Bx.gif")
        VizUtils.quickanim(By, timesteps, 10, vis_path_prefix * "By.gif")
        VizUtils.quickanim(u, timesteps, 10, vis_path_prefix * "u.gif")
        VizUtils.quickanim(v, timesteps, 10, vis_path_prefix * "v.gif")
    end

    return (x,y,u,v,Bx,By,p)
end

end