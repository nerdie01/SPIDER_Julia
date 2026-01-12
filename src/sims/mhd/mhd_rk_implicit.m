function [x,y,t,u,v,p,Bx,By] = mhd_rk_implicit(_omega, _A, _forcing, n, dt, timesteps, nu, eta, lorentz, skip, butcher, vis)
    if nargin < 12 || isempty(vis)
        vis = false;
    end

    addpath("../");

    max_iter_gl = 1000;
    tol_gl = 1e-6;

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

        s = length(butcher.b);
        for i_sub = 1:skip
            # initial guess
            [K_omega0, K_A0] = mhd_terms(omega, A, kx, ky, to_u, to_v, nu, eta, lorentz, forcing, mask);
            K_omega = repmat(K_omega0, [1, 1, s]);
            K_A = repmat(K_A0, [1, 1, s]);

            converged = false;
            for iter = 1:max_iter_gl
                omega_stages = zeros(n, n, s) + omega;
                A_stages = zeros(n, n, s) + A;

                for j = 1:s
                    for l = 1:s
                        omega_stages(:,:,j) = omega_stages(:,:,j) + dt * butcher.A(j,l) * K_omega(:,:,l);
                        A_stages(:,:,j) = A_stages(:,:,j) + dt * butcher.A(j,l) * K_A(:,:,l);
                    end
                end

                K_omega_new = zeros(n, n, s);
                K_A_new = zeros(n, n, s);
                for j = 1:s
                    [K_omega_new(:,:,j), K_A_new(:,:,j)] = mhd_terms(...
                        omega_stages(:,:,j), A_stages(:,:,j), ...
                        kx, ky, to_u, to_v, nu, eta, lorentz, forcing, mask);
                end

                err = 0;
                for j = 1:s
                    norm_omega = max(norm(K_omega(:,:,j), 'fro'), 1e-12);
                    norm_A = max(norm(K_A(:,:,j), 'fro'), 1e-12);
                    err = err + norm(K_omega_new(:,:,j) - K_omega(:,:,j), 'fro') / norm_omega;
                    err = err + norm(K_A_new(:,:,j) - K_A(:,:,j), 'fro') / norm_A;
                end
                err = err / (2*s);
                
                if err < tol_gl
                    K_omega = K_omega_new;
                    K_A = K_A_new;
                    converged = true;
                    break;
                end

                K_omega = K_omega_new;
                K_A = K_A_new;
            end
            
            if ~converged
                warning('fixed point iteration did not converge', t, i_sub, err);
            end

            omega_new = omega;
            A_new = A;
            for j = 1:s
                omega_new = omega_new + dt * butcher.b(j) * K_omega(:,:,j);
                A_new = A_new + dt * butcher.b(j) * K_A(:,:,j);
            end
            omega = omega_new;
            A = A_new;
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