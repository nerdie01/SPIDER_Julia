function [nl_omega, nl_A] = mhd_terms(omega, A, kx, ky, to_u, to_v, nu, eta, lorentz, forcing, mask)
    u = real(ifft2(to_u .* omega));
    v = real(ifft2(to_v .* omega));

    Bx = real(ifft2(1i*ky .* A));
    By = real(ifft2(-1i*kx .* A));

    wx = real(ifft2(1i*kx .* omega));
    wy = real(ifft2(1i*ky .* omega));
    
    adv_omega = fft2(-u.*wx - v.*wy) .* mask;

    curr = real(ifft2(-(kx.^2 + ky.^2) .* A));
    lorentz_omega = -lorentz * (1i*kx.*fft2(curr.*Bx) + 1i*ky.*fft2(curr.*By)) .* mask;

    visc_omega = -nu * (kx.^2 + ky.^2) .* omega;

    nl_omega = adv_omega + lorentz_omega + visc_omega + forcing;

    Ax = real(ifft2(1i*kx .* A));
    Ay = real(ifft2(1i*ky .* A));

    adv_A = fft2(-u.*Ax - v.*Ay) .* mask;
    res_A = -eta * (kx.^2 + ky.^2) .* A;

    nl_A = adv_A + res_A;
end