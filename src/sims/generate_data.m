clear;
pkg load symbolic;

addpath('mhd')
addpath('navier_stokes')


n = 128;
dt = 0.01;
timesteps = 100;
nu = 0.01;
eta = 0.02;
skip = 10;

butcher_rk4=[0, 0, 0, 0, 0; 0.5, 0.5, 0, 0, 0; 0.5, 0, 0.5, 0, 0; 1, 0, 0, 1, 0; 0, 1/6, 1/3, 1/3, 1/6];

omega = @(x, y) 3*sin(3*x).*sin(3*y) + sin(x).*sin(y);
A = @(x, y) cos(x) .* cos(y);
forcing = @(x, y) 0;

[x, y, t, u, v, p, Bx, By] = mhd_rk_generator(omega, A, forcing, n, dt, timesteps, nu, eta, 1, skip, butcher_rk4, true);
result = struct("x",x,"y",y,"t",t,"u",u,"v",v,"p",p,"Bx",Bx,"By",By)
save("-hdf5", "mhd_rk4_naive_result.hdf5", "result")