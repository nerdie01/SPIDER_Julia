module SpectralTerms

using FFTW

function mhd_terms(ω::AbstractMatrix, A::AbstractMatrix, ∂x::AbstractMatrix, ∂y::AbstractMatrix, uhat::AbstractMatrix, vhat::AbstractMatrix, ν::AbstractFloat, η::AbstractFloat, forcing::AbstractMatrix, mask::AbstractMatrix)::Tuple{Matrix{ComplexF64}, Matrix{ComplexF64}}
    ∇2 = ∂x.^2 + ∂y.^2
    curr = real(ifft(-∇2 .* A))

    u = real(ifft(uhat))
    v = real(ifft(vhat))
    Bx = real(ifft(∂y.*A))
    By = real(ifft(-∂x.*A))
    ωx = real(ifft(∂x.*ω))
    ωy = real(ifft(∂y.*ω))

    adv_ω = fft(-u.*ωx .- v.*ωy) .* mask
    lorentz_ω = (∂x.*fft(curr.*Bx) .+ ∂y.*fft(curr.*By)) .* mask
    visc_ω = ν.*∇2.*ω
    nl_ω = adv_ω .+ lorentz_ω .+ visc_ω .+ fft(forcing)

    Ax = real(ifft(∂x.*A))
    Ay = real(ifft(∂y.*A))
    adv_A = fft(-u.*Ax .- v.*Ay) .* mask
    res_A = η.*∇2.*A
    nl_A = adv_A .+ res_A

    return nl_ω, nl_A
end

end