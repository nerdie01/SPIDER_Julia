function [x,y,t,u,v,p,Bx,By] = mhd_rk_explicit(_omega, _A, _forcing, n, dt, timesteps, nu, eta, lorentz, skip, butcher, vis)
    if nargin < 12 || isempty(vis)
        vis = false
    end

    addpath("../")

    grid1d = (0:(n-1))/n*2*pi;
    [x, y] = meshgrid(grid1d);

    omega = _omega(x,y);
    A = _A(x,y);
    forcing = _forcing(x,y);

    u = zeros(n, n, timesteps);
    v = zeros(n, n, timesteps);
    p = zeros(n, n, timesteps);
    Bx = zeros(n, n, timesteps);
    By = zeros(n, n, timesteps);

    k = 0:n-1;
    k(k>n/2) = k(k>n/2) - n;

    mask = abs(k) <= n/3;
    mask = mask & mask.';

    kx = k;
    ky = k.';

    kk = kx.^2 + ky.^2;
    kk(1,1) = 1;

    to_u = 1i*ky ./ kk;
    to_v = -1i*kx ./ kk;
    to_p = 1 ./ kk;

    to_u(1,1) = 0;
    to_v(1,1) = 0;
    to_p(1,1) = 0;

    to_u = to_u .* mask;
    to_v = to_v .* mask;
    to_p = to_p .* mask;

    kx = kx .* mask;
    ky = ky .* mask;

    omega = fft2(omega) .* mask;
    A = fft2(A) .* mask;

    for t = 1:timesteps
        u(:,:,t) = real(ifft2(to_u .* omega));
        v(:,:,t) = real(ifft2(to_v .* omega));
        Bx(:,:,t) = real(ifft2(1i*ky .* A));
        By(:,:,t) = real(ifft2(-1i*kx .* A));

        ux = real(ifft2(1i*kx.*to_u.*omega));
        uy = real(ifft2(1i*ky.*to_u.*omega));
        vx = real(ifft2(1i*kx.*to_v.*omega));
        vy = real(ifft2(1i*ky.*to_v.*omega));

        Bxx = real(ifft2(-kx.*ky.*A));
        Bxy = real(ifft2(-ky.^2.*A));
        Byx = real(ifft2(kx.^2.*A));
        Byy = real(ifft2(ky.*kx.*A));

        p(:,:,t) = real(ifft2(to_p.*fft2(ux.^2 + vy.^2 + 2*uy.*vx - Bxx.^2 - Byy.^2 - 2*Bxy.*Byx) .* mask));

        f = @(x, y) mhd_terms(x, y, kx, ky, to_u, to_v, nu, eta, lorentz, forcing, mask);

        implicit = false;
        if norm(diag(butcher)) > 1e-12
            implicit = true;
        end

        for i = 1:skip
            [omega, A] = transform_fields_explicit(f, omega, A, dt, butcher);
        end

        if vis
            subplot(1,3,1);
            imagesc(real(ifft2(omega)));
            axis square;
            colorbar();
            title("Vorticity");

            subplot(1,3,2);
            imagesc(real(ifft2(A)));
            axis square;
            colorbar();
            title("B field potential");

            subplot(1,3,3);
            imagesc(p(:,:,t));
            axis square;
            colorbar();
            title("Pressure");

            drawnow;
        end
    end

    x = grid1d;
    y = grid1d;
    t = (0:timesteps-1)*dt*skip;
end

function [omega_P, A_P] = transform_fields_explicit(f, omega, A, dt, butcher)
    b_size = size(butcher, 1);

    k_omega = zeros(size(omega,1), size(omega,2), b_size - 1);
    k_A = zeros(size(A,1), size(A,2), b_size - 1);

    [k_omega(:,:,1), k_A(:,:,1)] = f(omega, A);
    k_omega(:,:,1) = dt * k_omega(1);
    k_A(:,:,1) = dt * k_A(1);

    a_omega_sum = zeros(size(omega));
    a_A_sum = zeros(size(A));

    d_omega = zeros(size(omega));
    d_A = zeros(size(A));

    for s = 2:b_size
        for i = 2:s
            a_omega_sum = a_omega_sum + butcher(s,i) * k_omega(:,:,s-1);
            a_A_sum = a_A_sum + butcher(s,i) * k_A(:,:,s-1);
        end

        [k_omega(:,:,s), k_A(:,:,s)] = f(omega + a_omega_sum, A + a_A_sum);
        k_omega(:,:,s) = dt * k_omega(:,:,s);
        k_A(:,:,s) = dt * k_A(:,:,s);

        d_omega = d_omega + butcher(b_size, s) * k_omega(:,:,s);
        d_A = d_A + butcher(b_size, s) * k_A(:,:,s);
    end

    omega_P = omega + d_omega;
    A_P = A + d_A;
end